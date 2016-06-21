'''
This file contains a series of scripts for evaluating gene models and determining their usability as
a "master gene"
'''

import pickle
from process_BC_data import gene_profile
from sklearn.mixture import GMM

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
    X = model.sample(NUM_TESTS)
    C = model.predict_proba(X)

    print(C)

    return "DONE!"

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
            print(calculate_shape_balance(gene_profiles[name]))