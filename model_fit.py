'''
Contains methods in order to pre-process and download data

Then to create and fit a gaussian model of
'''

import timeit
import pickle
import sqlite3
from process_BC_data import expression_profile
from process_BC_data import clinical_data

class sample:
    def __init__(self, profiles, sample_num):
        self.profiles = profiles
        self.sample_num = sample_num
'''
Creates sample containing all probes linked to a gene and dumps the array of samples into a pickle
'''
def get_expression_profiles(file):
    connection = sqlite3.connect("GeneExpression.db")
    cursor = connection.cursor()

    #first get a list of all the unique ids
    cursor.execute("Select Distinct Sample From BC_GeneExpression")
    rows = cursor.fetchall()
    sample_nums = [row[0] for row in rows]

    samples = []
    for num in sample_nums:
        #get all gene profiles from each profile
        cursor.execute("Select * From BC_GeneExpression WHERE Sample='%s' AND Gene != ''" % num)
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
def get_clinical_profiles(file):
    connection = sqlite3.connect("GeneExpression.db")
    cursor = connection.cursor()

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

if __name__ == "__main__":
    start = timeit.default_timer()

    print("preparting to dump pickles...")
    get_expression_profiles("BC_expression_profiles.pkl")
    get_clinical_profiles("BC_clinical_profiles.pkl")

    end = timeit.default_timer()

    print("DONE! Took " + str(end - start) + "s")
