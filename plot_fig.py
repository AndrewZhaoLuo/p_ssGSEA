'''
This file contains methods for plotting the evaluation graphs for each master_gene
'''

import os
import pickle
import sys
import cache_codec
import matplotlib

#Use backend which does not require X11 forwarding. Allows plotting on servers
matplotlib.use('Agg')

from matplotlib import pyplot as plotter

def get_ranking_data(enrichment_ranks):
    '''
    Given a series of enrichment rankings, returns a map of master gene names to a list with a number of elements
    equaling the total trials. Each element is a list with elements either a 1 or 0, 1 representing a gene set with the master gene and 0
    otherwise

    :param enrichment_ranks: the path and file containing expression information to read
    :type enrichment_ranks: str

    :return: a mapping of master genes to a list or lists
    '''

    gene_sets = cache_codec.load_filtered_gene_sets('BC')

    ranks = {}

    for gene in enrichment_ranks:
        print(gene)
        ranks[gene] = []
        p_values_gene = enrichment_ranks[gene]

        good_sets = {1}
        good_sets.clear()
        for set in gene_sets:
            if gene in gene_sets[set].genes:
                good_sets.add(set)

        for trial in p_values_gene:
            trial_rank = []

            #iterator that goes through sets in increasing order of tstat
            sorted_gene_sets = sorted(trial, key=trial.get, reverse=True)
            for set in sorted_gene_sets:
                if set in good_sets:
                   trial_rank.append(1)
                else:
                    trial_rank.append(0)

            ranks[gene].append(trial_rank)

    return ranks

def make_graph_against_sets_picked_top(gene, rankings, x_range, method):
    '''
    Graphing overall performance across all master_genes
    '''
    gene_rank = rankings[gene]
    pop = sum(gene_rank[0])
    y_axis = []
    for x in x_range:
        y_local = 0
        y_tot = 0
        for trial in gene_rank:
            picked_sets = sum(trial[:(x + 1)])

            if picked_sets > 0:
                y_local += 1
            y_tot += 1
        y_axis.append(y_local / y_tot)

    x_axis = [x + 1 for x in x_range]
    plotter.axis([0, len(x_axis), 0, 1.0])
    plotter.xlabel("Top # of sets picked")
    plotter.ylabel("Proportion of Trials With True Gene Sets")
    plotter.title(gene)

    plotter.plot(x_axis, y_axis)
    plotter.savefig("./Pictures/" + method + "Analysis_" + gene + " effectiveness: " + str(pop))
    plotter.close()

def make_graph_against_sets_all_top(rankings, x_range, method):
    '''
    Graphing performance for a single master gene
    '''
    y_axis = []
    for x in x_range:
        y_tot_hits = 0
        y_tot_sets = 0

        for gene in rankings:
            gene_rank = rankings[gene]
            for trial in gene_rank:
                picked_sets = sum(trial[:(x + 1)])

                if picked_sets > 0:
                    y_tot_hits += 1
                y_tot_sets += 1

        y_axis.append(y_tot_hits / y_tot_sets)

    print(y_axis)
    x_axis = [x + 1 for x in x_range]
    plotter.plot(x_axis, y_axis, label='mean')

    plotter.axis([0, len(x_axis), 0, 1])
    plotter.xlabel("Top # of sets picked")
    plotter.ylabel("Proportion of True Gene Sets")
    plotter.savefig("./Pictures/" + str(method) + "TotalAnalysis")
    plotter.close()


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "":
        print("Give a test set result!")
        exit()

    test = sys.argv[1]

    #This loading just happens to be where we dump the pickle.
    enrichment_ranks = pickle.load(open(os.getcwd() + '/Data/AppCache/BC/' + test +  'CachedEnrichmentPValueSplit.pkl', 'rb'))

    ranking = get_ranking_data(enrichment_ranks)

    all_genes = ranking.keys()
    print(len(all_genes))
    for gene in all_genes:
        make_graph_against_sets_picked_top(gene, ranking, range(0, 100), test)

    make_graph_against_sets_all_top(ranking, range(0, 50), test)
