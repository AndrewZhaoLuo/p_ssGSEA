'''
This file contains methods for plotting various graphs given various data points
'''

import os
import pickle
import cache_codec

enrichment_ranks = pickle.load(open(os.getcwd() + '/Data/AppCache/BC/CachedEnrichmentPValueSplit.pkl', 'rb'))

#list of lists. each element is a list of 0's and 1's. 1 means the master gene was in the set at rank indicie + 1
def get_ranking_data(enrichment_ranks):
    gene_sets = cache_codec.load_all_gene_sets()

    ranks = []

    for gene in enrichment_ranks:
        p_values_gene = enrichment_ranks[gene]
        for trial in p_values_gene:
            trial_rank = []

            #iterator that goes through sets in increasing order of tstat
            sorted_gene_sets = sorted(trial, key=trial.get)
            for set in sorted_gene_sets:
                if gene in gene_sets[set].genes:
                   trial_rank.append(1)
                else:
                    trial_rank.append(0)

            ranks.append(trial_rank)

    return ranks