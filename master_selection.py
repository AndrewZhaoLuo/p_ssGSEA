'''
This file contains a series of scripts for evaluating gene models and determining their usability as
a "master gene"
'''

import timeit
import pickle
import sqlite3
from process_BC_data import gene_set
from process_BC_data import gene_profile
from random import gauss

import heapq

'''
returns the π measure as per the paper
'''
def calculate_prior(model):
    coeffs = model.weights_

    return min(coeffs[0], 1- coeffs[0])
'''
Todo: use an analytical method

Unsure variance in sampling method
'''
def calculate_bayes_error(model):
    coeffs = model.weights_
    mus = [x[0] for x in model.means_]
    sigmas = [x[0] ** 0.5 for x in model.covars_]

    '''
    #old sampling method
    NUM_TESTS = 10000
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
    '''

def calculate_fold_change(model):
    mus = model.means_
    return abs(mus[0] - mus[1])

def calculate_shape_balance(model):
    sigmas = model.covars_

    return max(sigmas[1] / sigmas[0], sigmas[0] / sigmas[1])

def calculate_popularity(name):
    connection = sqlite3.connect("GeneExpression.db")
    cursor = connection.cursor()
    cursor.execute("Select GeneSet From BC_GeneSet_Genes WHERE Gene='%s'" % name)

    return len(cursor.fetchall())

'''
Model is valid if prior >= 0.1, non-zero popularity.
Then divide priors of [0.1, 0.5] into 10 bins for each bin take the best (lowest bayes error) 10
genes for each bin and dump
'''
def dump_best_models(gene_profiles, num_bins, genes_per_bin):
    bin_bounds = [0.1] + [0.1 + x * (0.5 - 0.1) / num_bins for x in range(1, num_bins + 1)]
    find_bin = lambda prior: num_bins if prior == 0.5 \
                                else sum([i for i in range(0, len(bin_bounds) - 1)
                                  if (prior >= bin_bounds[i] and prior < bin_bounds[i + 1])])

    bins = [{} for x in range(0, num_bins)]

    start = timeit.default_timer()

    #Here we check every gene fits the criteria and then add it to appropriate bin
    #bins are maps k : v where k = gene name v = bayes error
    i = 0
    gene_names = gene_profiles.keys()
    for gene in gene_names:
        model = gene_profiles[gene]

        popularity = calculate_popularity(gene)
        prior = calculate_prior(model)
        error = calculate_bayes_error(model)

        #pick and get all gene errors
        if popularity > 0 and prior >= 0.1:
            bin = find_bin(prior)
            bins[bin][gene] = error

        i += 1
        if i % 100 == 0:
            print("Completed classifying " + str(i) + " genes\t" + str(timeit.default_timer() - start) +"s")
            start = timeit.default_timer()

    #then we go through each bin and take the 10 smallest bayes error
    #bestmodels is k: v where k = gene name and v = bayes error for that model
    best_models = {}
    for bin in bins:
        best_genes = heapq.nsmallest(10, genes_per_bin, key=bin.get)
        print(bin)
        print(best_genes)
        print()
        for gene in best_genes:
            best_models[gene] = bin[gene]

    print(best_models)
    pickle.dump(best_models, open("BC_master_genes.pkl", 'wb'))


if __name__ == "__main__":
    master_genes = pickle.load(open("BC_master_genes_old.pkl", 'rb'))

    best_genes = sorted(master_genes, key=master_genes.get)
    for gene in best_genes:
        print(gene)
        print("\t" + str(master_genes[gene]))
        print()

    '''
    #code for loading master genes
    print("Loading data set")
    #gene_sets = pickle.load(open("BC_gene_sets.pkl",'rb'))
    gene_profiles = pickle.load(open("BC_trained_models.pkl", 'rb'))

    print("Examining genes...")

    start = timeit.default_timer()
    dump_best_models(gene_profiles)
    end = timeit.default_timer()

    print("DONE! Took " + str(end - start) + "s")
    '''