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
def dump_expression_profiles(file, query, cursor):
    #first get a list of all the unique ids
    cursor.execute("Select Distinct Sample From BC_GeneExpression")
    rows = cursor.fetchall()
    sample_nums = [row[0] for row in rows]

    samples = []
    for num in sample_nums:
        #get all gene profiles from each profile
        #cursor.execute("Select * From BC_GeneExpression WHERE Sample='%s' AND Gene != ''" % num)
        cursor.execute(query)
        sample_profiles = cursor.fetchall()

        profiles = []
        for profile in sample_profiles:
            sample_num, substance, gene, log_ratio, log_error, p_value, intensity, flag = profile

            profiles.append(expression_profile(sample_num, substance, gene, log_ratio,
                                               log_error, p_value, intensity, flag))
        samples.append(sample(profiles, num))

    pickle.dump(samples, open(file, 'wb'))

'''
Downloads clinical profiles and dumps the array of them into a pickle
'''
def dump_clinical_profiles(file, query, cursor):
    #cursor.execute("Select * From BC_ClinicalData")
    cursor.execute(query)
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

def read_dumped_data(file):
    return pickle.load(open(file, 'rb'))

if __name__ == "__main__":
    connection = sqlite3.connect("GeneExpression.db")
    cursor = connection.cursor()



    connection.commit()

