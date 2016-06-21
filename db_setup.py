'''
This is a simple series of script which

IN GENERAL THIS MESSES WITH THE DB SCHEMA, DONT MESS WITH UNLESS YOU KNOW WHAT YOU ARE DOING!!!
'''

import sqlite3
import process_BC_data as BC
import pickle

from process_BC_data import expression_profile
from process_BC_data import clinical_data
from process_BC_data import sample
from process_BC_data import gene_profile

import timeit

'''
DATABASE CODE
'''
def create_BC_schema(cursor):
    cursor.execute('''CREATE TABLE BC_GeneExpression
        (Sample int, Substance text, Gene text, Log_Ratio real,
         Log_Ratio_Error real, p_Value real, Intensity real, Flag int)''')

    cursor.execute('''CREATE TABLE BC_ClinicalData
        (Sample int, FirstSeriesID int, Posnodes text, EventMeta int,
        EventDeath int, TimeSurvival real, TimeRecur real, TimeMeta real, ESR1 int,
        NIH int, StGallen int, Conserv int, C1FromData real, C1CrossValid real,
        C1Used real)''')

def load_BC_data(cursor):
    gene_profiles = BC.getExpressionProfiles(BC.BC_EXPRESSION_DIR)
    for profile in gene_profiles:
        cursor.execute("INSERT INTO BC_GeneExpression VALUES(?,?,?,?,?,?,?,?)",
                       (profile.sample_num, profile.substance, profile.gene,
                        profile.log_ratio, profile.log_error, profile.p_value,
                        profile.intensity, profile.flag))

    clinical_profiles = BC.getClinicalData(BC.BC_CLINICAL_DATA_FILE)
    for profile in clinical_profiles:
        cursor.execute("INSERT INTO BC_ClinicalData VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                       (profile.sample_num, profile.first_series, profile.posnodes, profile.event_meta,
                        profile.event_death, profile.time_survival, profile.time_recur, profile.time_meta,
                        profile.esr1, profile.nih, profile.st_gallen, profile.conserv,
                        profile.c1_from_data, profile.c1_cross_valid, profile.c1_used))

'''
Tears down and rebuild BC database (PLEASE DON'T CALL UNLESS YOU KNOW WHAT YOU ARE DOING!!!!!)
'''
def rebuild_BC_db(cursor):
    print("Rebuilding BC database...")
    start = timeit.default_timer()

    cursor.execute("DROP TABLE IF EXISTS BC_GeneExpression")
    cursor.execute("DROP TABLE IF EXISTS BC_ClinicalData")

    create_BC_schema(cursor)
    load_BC_data(cursor)

    end = timeit.default_timer()
    print("Finished building DB! Took " + str(end - start) + "s")


'''
GETTING DATA FROM DATABASE CODE
'''
'''
Creates sample containing all probes linked to a gene and dumps the array of samples into a pickle
'''
def dump_expression_profiles(file, cursor):
    #first get a list of all the unique ids
    cursor.execute("Select Distinct Sample From BC_GeneExpression")
    rows = cursor.fetchall()
    sample_nums = [row[0] for row in rows]

    samples = []
    for num in sample_nums:
        #get all gene profiles from each profile
        cursor.execute("Select * From BC_GeneExpression WHERE Sample='%s' AND Gene != ''" % num)
        #cursor.execute(query)
        sample_profiles = cursor.fetchall()

        profiles = []
        for profile in sample_profiles:
            sample_num, substance, gene, log_ratio, log_error, p_value, intensity, flag = profile

            profiles.append(expression_profile(sample_num, substance, gene, log_ratio,
                                               log_error, p_value, intensity, flag))
        samples.append(sample(profiles, num))

    pickle.dump(samples, open(file, 'wb'))

'''
Dump gene containing all probes linked to a gene and dumps the array of samples into a pickle
'''
def dump_gene_profiles(file, cursor):
    #first get a list of all the unique ids
    cursor.execute("Select Distinct Gene From BC_GeneExpression")
    rows = cursor.fetchall()
    gene_names = [row[0] for row in rows]

    genes = []
    for gene_name in gene_names:
        #get all gene intensities for each gene
        cursor.execute("Select Intensity From BC_GeneExpression WHERE Gene= '%s'" % gene_name)
        intensities_rows = cursor.fetchall()

        #importing negatives seems to append
        #intensities = [row[0] for row in intensities_rows]
        #intensities = [row[0] for row in intensities_rows]
        print(gene_name + " " + str(intensities_rows))
        #genes.append(gene(gene_name, intensities))

    pickle.dump(genes, open(file, 'wb'))

'''
Downloads clinical profiles and dumps the array of them into a pickle
'''
def dump_clinical_profiles(file, cursor):
    cursor.execute("Select * From BC_ClinicalData")
    rows = cursor.fetchall()

    profiles = []
    for row in rows:
        sample_num, first_series, posnodes, event_meta, \
        event_death, time_survival, time_recur, time_meta, esr1, nih, st_gallen, conserv, \
        c1_from_data, c1_cross_valid, c1_used = row
        profiles.append(clinical_data(sample_num, first_series, posnodes, event_meta, \
                                      event_death, time_survival, time_recur, time_meta,
                                      esr1, nih, st_gallen, conserv, \
                                      c1_from_data, c1_cross_valid, c1_used))

    pickle.dump(profiles, open(file, 'wb'))

'''
Given a pickle of samples, dumps into gene_sample form
'''
def filter_sample_profile_to_gene(file):
    samples = []
    for profile in read_dumped_data("BC_expression_profiles.pkl"):
        samples.append(profile)

    genes = []
    num_genes = len(samples[0].profiles)
    for i in range(0, num_genes):
        gene_name = samples[0].profiles[i].gene
        intensities = []
        sample_nums = []
        for person in samples:
            intensities.append(person.profiles[i].intensity)
            sample_nums.append(person.sample_num)
            #if gene_name == "ERBB2":
            #    print(person.profiles[i].gene)
            #    print("\t" + str(person.profiles[i].intensity))
            #print([x.gene for x in person.profiles])
        genes.append(gene_profile(intensities, sample_nums, gene_name))

    check_length = len(genes[0].intensities)

    #check this method works
    for gene in genes:
        assert check_length == len(gene.intensities)

    pickle.dump(genes, open(file, 'wb'))

def read_dumped_data(file):
    return pickle.load(open(file, 'rb'))

if __name__ == "__main__":
    connection = sqlite3.connect("GeneExpression.db")
    cursor = connection.cursor()

    filter_sample_profile_to_gene("BC_gene_profiles.pkl")

    connection.commit()

