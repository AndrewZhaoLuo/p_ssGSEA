import pickle

from master_selection import calculate_popularity
from data_models import *

import sqlite3

#ToDo: documentation!

'''
Downloads clinical data files and dumps to pickle
'''
def dump_clinical_profiles(file, cursor):
    cursor.execute("Select * From BC_ClinicalData")
    rows = cursor.fetchall()

    profiles = []
    for row in rows:
        sample_num, first_series, posnodes, event_meta, \
        event_death, time_survival, time_recur, time_meta, esr1, nih, st_gallen, conserv, \
        c1_from_data, c1_cross_valid, c1_used = row

        other_fields = {}
        other_fields["first_series"] = first_series
        other_fields["posnodes"] = posnodes
        other_fields["event_meta"] = event_meta
        other_fields["event_death"] = event_death
        other_fields["time_survival"] = time_survival
        other_fields["time_recur"] = time_recur
        other_fields["time_meta"] = time_meta
        other_fields["esr1"] = esr1
        other_fields["nih"] = nih
        other_fields["st_gallen"] = st_gallen
        other_fields["conserv"] = conserv
        other_fields["c1_from_data"] = c1_from_data
        other_fields["c1_cross_valid"] = c1_cross_valid
        other_fields["c1_used"] = c1_used

        profiles.append(clinical_data(sample_num, other_fields))

    pickle.dump(profiles, open(file, 'wb'))

'''
Downloads expression data files and dumps to pickle
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

        profiles = {}
        for profile in sample_profiles:
            sample_num, substance, gene, log_ratio, log_error, p_value, intensity, flag = profile
            other_fields = {}
            other_fields["log_ratio"] = log_ratio
            other_fields["log_error"] = log_error
            other_fields["p_value"] = p_value
            other_fields["flag"] = flag
            other_fields["substance"] = substance

            profiles[gene] = (expression_profile(sample_num, gene, intensity, other_fields))

        samples.append(sample(profiles, num))

    pickle.dump(samples, open(file, 'wb'))

'''
precompute the popularity of all genes for future use
'''
def dump_gene_popularity(gene_names):
    gene_pop = {}

    for names in gene_names:
        gene_pop[names] = calculate_popularity(names)

    pickle.dump(gene_pop, open("BC_gene_popularity.pkl", 'wb'))

'''
dumps all samples man, maps sample_num to sample objects
'''
def dump_gene_samples(file, cursor):
    #first get all unique id's from the db
    cursor.execute("Select Distinct Sample From BC_ClinicalData")
    rows = cursor.fetchall()
    ids = [row[0] for row in rows]

    samples = {}
    #for each - create the sample
    for id in ids:
        cursor.execute("Select * From BC_GeneExpression WHERE Sample='%s' AND NOT Gene=''" % id)
        rows = cursor.fetchall()

        profiles = {}
        for row in rows:
            sample_num, substance, gene, log_ratio, log_error, p_value, intensity, flag = row

            other_fields = {}
            other_fields["log_ratio"] = log_ratio
            other_fields["log_error"] = log_error
            other_fields["p_value"] = p_value
            other_fields["flag"] = flag
            other_fields["substance"] = substance

            profile = expression_profile(sample_num, gene, intensity, other_fields)
            profiles[gene] = profile


        samples[id] = sample(profiles, id)
        print('yay!')

    pickle.dump(samples, open(file, 'wb'))