'''
This file contains script for emulating the ssGSEA algorithm from Barbie et al. 2009
'''

import model_fit
import numpy as np
from simulation import *

def calculate_enrichment_score(gene_set, expressions, omega):
    """
    Given a gene set, a map of gene names to expression levels, and a weight omega, returns the ssGSEA
    enrichment score for the gene set as described by *D. Barbie et al 2009*

    :requires: every member of gene_set is a key in expressions

    :param gene_set: a set of gene_names in the set
    :type gene_set: set

    :param expressions: a dictionary mapping gene names to their expression values
    :type expressions: dict

    :param omega: the weighted exponent on the :math:`P^W_G` term.
    :type omega: float

    :returns: an array representing the intermediate Enrichment Scores for each step along the sorted gene list.
              To find the total enrichment score, take the sum of all values in the array. To replicate a result
              result closer to GSEA, take the max of the values in the array
    """
    #list of keys in order to parse to get low -> high values
    keys_sorted = sorted(expressions, key=expressions.get)

    #running summations of the PGW and PNG scores in the Barbie et al 2009 paper's algo. desc.
    P_GW_numerator = 0
    #dead_genes = [x for x in gene_set if x not in expressions.keys()]
    #print(dead_genes)
    P_GW_denominator = sum([abs(expressions[gene]) ** omega for gene in gene_set])
    P_GW = lambda : P_GW_numerator / P_GW_denominator

    P_NG_numerator = 1
    P_NG_denominator = len(expressions) - len(gene_set)
    P_NG = lambda : P_NG_numerator / P_NG_denominator

    scores = []
    for gene in keys_sorted:
        #hitting a gene in the gene set
        if gene in gene_set:
            P_GW_numerator += abs(expressions[gene]) ** omega
        else:
            P_NG_numerator += 1

        scores.append(P_GW() - P_NG())

    return scores

#Below methods are meant to download and dump scores and shit man

'''
ToDo: migrate and cleanup, dumps pickles with models based on the enrichments scores of each data set
'''
def get_model_pathways():
    """
    ToDo: Move to caching
    """
    gene_sets = pickle.load(open("BC_good_gene_sets.pkl", 'rb'))
    samples = pickle.load(open("BC_sample_profiles.pkl", 'rb'))
    print("Data loaded!")

    models = {}
    #for each gene set
    for set in gene_sets:
        gene_set = set.genes

        #go through all the samples and calculate the ES
        scores = []
        for id in samples.keys():
            profile = samples[id].profiles
            expressions = {}

            for gene in profile.keys():
                expressions[gene] = profile[gene].intensity

            score = calculate_enrichment_score(gene_set, expressions, 0.25)
            scores.append(sum(score))

        #add the models bro!
        model = model_fit.fit_test_model([[x] for x in scores])
        models[set.set_name] = model
        print("Set " + str(set.set_name) + " done!")

    pickle.dump(models, open("BC_enriched_set_models.pkl", 'wb'))

'''
ToDo: make extensible and reusable

simulates phenotype assuming given gene is master gene.
Then plots the gaussian distribution of the enrichment scores for the given gene set
broken up by phenotype

#on gene sets
'''
def ssgsea_on_simulation(gene, models, profiles, gene_set, title):
    """
    ToDo: break up into methods to get the scores, and another to plot
    """

    #simulate phenotype labels and seperate the class0's/class1's
    model = models[gene]
    phenotypes = simulate_data(model, gene, profiles, len(profiles))

    #these contain the expressions for those in class0/class1 respectively
    class0 = []
    class1 = []
    for profile in profiles:
        expression_p = profile.profiles
        expressions = {}
        for key in expression_p:
            expressions[key] = expression_p[key].intensity
        expression = sum(calculate_enrichment_score(gene_set, expressions, 0.25)) / len(gene_set)
        key = profile.sample_num
        if phenotypes[key] == 1:
            class1.append(expression)
        else:
            class0.append(expression)

    #we then want the expression profiles of the ERBB2 gene to work correctly
    #creating models!
    sigma0 = np.std(class0)
    mu0 = np.mean(class0)

    sigma1 = np.std(class1)
    mu1 = np.mean(class1)

    #gaussian_sampling.plot_multidist([10000, 10000], [mu0, mu1], [sigma0, sigma1], [1, 1], "1_" + title + "_cleaned", False)
    #gaussian_sampling.plot_multidist_from_values(title="2_" + title + "_raw", values=[class0, class1])
    score = abs(mu0 - mu1) / (sigma0 + sigma1)
    print('tic!')
    return score

if __name__ == "__main__":
    #first load all the models you want
    gene_sets = pickle.load(open("BC_good_gene_sets.pkl", 'rb'))
    models = pickle.load(open("BC_trained_models.pkl", 'rb'))
    profiles = pickle.load(open("BC_expression_profiles.pkl", 'rb'))

    print("Loaded models!")

    for set in gene_sets:
        ssgsea_on_simulation("ERBB2", models, profiles, set.genes, set.set_name)




