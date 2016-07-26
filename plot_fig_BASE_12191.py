'''
This file contains methods for plotting various graphs given various data points
'''

import os
import pickle
import cache_codec
import numpy as np
from matplotlib import pyplot as plotter
matplotlib.use('Agg')

enrichment_ranks = pickle.load(open(os.getcwd() + '/Data/AppCache/BC/ssGSEACachedEnrichmentPValueSplit.pkl', 'rb'))

#list of lists. each element is a list of 0's and 1's. 1 means the master gene was in the set at rank indicie + 1
def get_ranking_data(enrichment_ranks):
    gene_sets = cache_codec.load_all_gene_sets()

    ranks = {}

    for gene in enrichment_ranks:
        ranks[gene] = []
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

            ranks[gene].append(trial_rank)

    return ranks

def make_graph_against_sets_picked(gene, rankings, x_range):
    '''
    currently averages all sets for a single gene
    '''

    gene_rank = rankings[gene]
    pop = sum(gene_rank[0])
    y_best = []
    y_worst = []
    y_median = []
    y_mean = []
    for x in x_range:
        y_local = []
        for trial in gene_rank:
            tot_sets = sum(trial)
            picked_sets = sum(trial[:(x + 1)])
            percent_picked = picked_sets / tot_sets
            y_local.append(percent_picked)

        y_mean.append(sum(y_local) / len(y_local))
        y_median.append(np.median(y_local))
        y_best.append(max(y_local))
        y_worst.append(min(y_local))

    x_axis = [x for x in x_range]
    plotter.xlabel("Top # of sets picked")
    plotter.ylabel("Proportion of True Gene Sets")
    plotter.title(gene)
    plotter.plot(x_axis, y_mean, label='mean')
    plotter.plot(x_axis, y_median, label='median')
    plotter.plot(x_axis, y_best, label='best')
    plotter.plot(x_axis, y_worst, label='worse')
    plotter.legend(loc='upper left')

    plotter.savefig("./Pictures/Analysis_" + gene + " popularity: " + str(pop))
    plotter.close()

def make_graph_against_sets_all(rankings, x_range):
    '''
    currently averages all sets across all genes
    '''
    pop = 0

    y_tot_sets = []
    y_tot_hit = []
    for x in x_range:
        y_tot_sets_local = []
        y_tot_hit_local = []

        for gene in rankings:
            gene_rank = rankings[gene]
            for trial in gene_rank:
                tot_sets = sum(trial)
                picked_sets = sum(trial[:(x + 1)])

                y_tot_sets_local.append(tot_sets)
                y_tot_hit_local.append(picked_sets)

        y_tot_sets.append(y_tot_sets_local)
        y_tot_hit.append(y_tot_hit_local)

    y_mean = []
    for i in range(0, len(y_tot_hit)):
        y_mean.append(sum(y_tot_hit[i]) / sum(y_tot_sets[i]))

    x_axis = [x for x in x_range]
    plotter.plot(x_axis, y_mean, label='mean')

    plotter.xlabel("Top # of sets picked")
    plotter.ylabel("Proportion of True Gene Sets")
    plotter.savefig("./Pictures/TotalAnalysis")
    plotter.close()


ranking = get_ranking_data(enrichment_ranks)
make_graph_against_sets_picked("ERBB2", ranking, range(0, len(ranking["ERBB2"][0])))
make_graph_against_sets_picked("SOX10", ranking, range(0, len(ranking["SOX10"][0])))
make_graph_against_sets_all(ranking, range(0, len(ranking["SOX10"][0])))
