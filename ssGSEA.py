'''
This file contains script for emulating the ssGSEA algorithm from Barbie et al. 2009
'''

import pickle
import matplotlib.pyplot as plotter
from data_models import gene_set

'''
Calculates the enrichment score of this profile

gene_set = list of genes in this set, every gene in the gene set should be in samples as a key value

samples = a mapping of gene names -> intensity values

omega = the weight bro of the P_G^W term. 0 = just like that russian stat test

returns an array representing the intermediate ES values for each step along the random walk
'''
def calculate_enrichment_score(gene_set, expressions, omega):
    #list of keys in order to parse to get low -> high values
    keys_sorted = sorted(expressions, key=expressions.get)

    #running summations of the PGW and PNG scores in the Barbie et al 2009 paper's algo. desc.
    P_GW_numerator = 0
    P_GW_denominator = sum([abs(expressions[gene]) ** omega for gene in gene_set])
    P_GW = lambda : P_GW_numerator / P_GW_denominator

    P_NG_numerator = 1
    P_NG_denominator = len(expressions) - len(gene_set)
    P_NG = lambda : P_NG_numerator / P_NG_denominator

    #print("P_GW\tP_NG")
    scores = []
    ES = 0
    #i = 0
    for gene in keys_sorted:
        #hitting a gene in the gene set
        if gene in gene_set:
            P_GW_numerator += abs(expressions[gene]) ** omega
        else:
            P_NG_numerator += 1
        #i += 1

        #print(str(P_GW()) + "\t\t" + str(P_NG()))

        ES += P_GW() - P_NG()
        scores.append(P_GW() - P_NG())

    return scores

if __name__ == "__main__":
    gene_set = pickle.load(open("BC_gene_sets.pkl", 'rb'))
    samples = pickle.load(open("BC_sample_profiles.pkl", 'rb'))
    print("loaded!")

    profile = samples[4].profiles
    gene_set = [gene for gene in gene_set[2].genes if gene in profile.keys()]
    print(gene_set)

    expressions = {}
    for gene in profile.keys():
        expressions[gene] = profile[gene].intensity

    scores = calculate_enrichment_score(gene_set, expressions, 2)
    plotter.plot([i + 1 for i in range(0, len(scores))], scores)
    plotter.savefig("SSGSEA.png")