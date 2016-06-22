'''
This file contains a series of scripts for evaluating gene models and determining their usability as
a "master gene"
'''

import pickle
from process_BC_data import gene_profile
from random import gauss


'''
returns the π measure as per the paper
'''
def calculate_prior(model):
    coeffs = model.weights_

    return min(coeffs[0], 1- coeffs[0])

'''
ToDo
'''
def calculate_bayes_error(model):
    NUM_TESTS = 1000000
    mus = [x[0] for x in model.means_]
    sigmas = [x[0] ** 0.5 for x in model.covars_]

    #print(mus)
    #print(sigmas)
    #actual values from the paper
    #mus = [1.08, 1.65]
    #sigmas = [0.14, 0.09]

    missClass0 = 0
    missClass1 = 0

    for i in range(0, NUM_TESTS):
        sample0 = gauss(mus[0], sigmas[0])
        zscore0 = abs(sample0 - mus[0]) / sigmas[0]
        zscore1 = abs(sample0 - mus[1]) / sigmas[1]
        missClass1 += int(zscore1 < zscore0)

        sample1 = gauss(mus[1], sigmas[1])
        zscore0 = abs(sample1 - mus[0]) / sigmas[0]
        zscore1 = abs(sample1 - mus[1]) / sigmas[1]
        missClass1 += int(zscore0 < zscore1)

    return (abs(missClass0) + abs(missClass1)) / (NUM_TESTS * 2)

def calculate_fold_change(model):
    mus = model.means_
    return abs(mus[0] - mus[1])

def calculate_shape_balance(model):
    sigmas = model.covars_

    return max(sigmas[1] / sigmas[0], sigmas[0] / sigmas[1])

if __name__ == "__main__":
    gene_profiles = pickle.load(open("BC_trained_models.pkl", 'rb'))

    gene_names = gene_profiles.keys()
    for name in gene_names:
        if name == "ERBB2":
            print(calculate_bayes_error(gene_profiles[name]))