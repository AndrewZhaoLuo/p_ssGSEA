'''
This file contains scripts for running stuff man
'''

import cache_codec
import numpy

from analyze_enrichment import rank_by_t_test
from analyze_enrichment import evaluate_rankings

def run_analysis_on_dataset(data_set, n):
    '''
    Given data_set, runs analysis on all genes as masters in ranks one which occur in multiple gene_sets
    '''
    master_genes = cache_codec.load_gene_popularity(data_set).keys()
    #for master_gene in master_genes:

    gene_evaluation = {}
    gene_evaluation_full = {}
    for master_gene in master_genes:
        if "/" not in master_gene and master_gene=="SOX10":
            print("Starting analysis on gene " + str(master_gene))
            enrichment_scores = cache_codec.load_ssGSEA_scores(data_set)
            phenotypes = cache_codec.load_sim_phenotypes(data_set, n, master_gene)

            rankings = rank_by_t_test(enrichment_scores, phenotypes)

            print(rankings)

            gene_sets = cache_codec.load_filtered_gene_sets(data_set)
            evaluation = evaluate_rankings(rankings, gene_sets, master_gene)

            if len(evaluation) > 0:
                print(master_gene + "\t" + str(evaluation))
                gene_evaluation[master_gene] = numpy.median(evaluation)
                gene_evaluation_full[master_gene] = evaluation

    sorted_genes = sorted(gene_evaluation, key=gene_evaluation.get)
    f = open('results.txt','w')
    f.write("Gene\tMedian Rank\tFull Ranks\n")
    for gene in sorted_genes:
        f.write(str(gene) + "\t" + str(gene_evaluation[gene]) + "\t" + str(gene_evaluation_full[gene]) + "\n")

    f.close()

if __name__ == "__main__":
    print("Starting analysis...")
    #redirect stdout to nothing in order to speed up program
    import os
    import sys

    f = open(os.devnull, 'w')
    sys.stdout = f

    run_analysis_on_dataset("BC", 10)