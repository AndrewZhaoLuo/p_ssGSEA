'''
This file contains scripts for running stuff man
'''

from multiprocessing import Pool

import numpy

#import cache_codec
from cache_codec import load_sim_phenotype_keyed
from cache_codec import load_gene_popularity
from cache_codec import load_ssGSEA_scores
from cache_codec import load_filtered_gene_sets

from analyze_enrichment import rank_by_t_test_keyed
from analyze_enrichment import evaluate_rankings_keyed

NUM_PROCESSES = 4

def run_analysis_on_dataset(data_set, n):
    '''
    Given data_set, and n which is the number of replicates for the experiment:
    for every gene in the data_set for which there exists a filtered gene set, creates imaginary phenotypes for every
    sample using a gaussian mixture model. Imaginary phenotypes are then used to fit gaussian models of enrichment values
    for each gene set. A gene set is considered "intersting" if the generated models can seperate. This is done
    by a t-test of the class sepearated enrichment data. Gene ranking is done by looking at the sets in which the
    master gene belongs. The master gene should show up at the very top of the sets.

    ToDo: finish
    '''

    master_genes = load_gene_popularity(data_set).keys()
    enrichment_scores = load_ssGSEA_scores(data_set)
    gene_sets = load_filtered_gene_sets(data_set)

    gene_evaluation_median = {}
    gene_evaluation_full = {}

    good_genes = []
    for master_gene in master_genes:
        #ToDo: pre-entry sanitation of the database!
        if "/" not in master_gene:
            good_genes.append(master_gene)

    #fast and threaded!!!
    pool = Pool(processes=NUM_PROCESSES)

    #load phenotypes, maps master gene to phenotype
    pheno_params = [(data_set, n, good_gene) for good_gene in good_genes]
    loaded_phenos = pool.starmap(load_sim_phenotype_keyed, pheno_params)

    #process into more usable form of master gene -> phenotypes
    pheno_map = {}
    for map in loaded_phenos:
        for key in map.keys():
            pheno_map[key] = map[key]

    #work out the enrichment_score ranks
    enrichment_list = pool.starmap(rank_by_t_test_keyed, [(enrichment_scores, pheno_map[gene], gene) for gene in good_genes])
    enrichment_map = {}
    for enrichment in enrichment_list:
        for key in enrichment.keys():
            enrichment_map[key] = enrichment[key]

    #finally translate that to gene rankings
    results_map = pool.starmap(evaluate_rankings_keyed, [(enrichment_map[gene], gene_sets, gene) for gene in good_genes])
    for result in results_map:
        for master_gene in result.keys():
            evaluation = result[master_gene]
            if len(evaluation) > 0:
                print(master_gene + "\t" + str(evaluation))
                gene_evaluation_median[master_gene] = numpy.median(evaluation)
                gene_evaluation_full[master_gene] = evaluation

    pool.close()
    pool.join()

    #print out final rankings
    sorted_genes = sorted(gene_evaluation_median, key=gene_evaluation_median.get)
    f = open('results.txt','w')
    f.write("Gene\tMedian Rank\tFull Ranks\n")
    for gene in sorted_genes:
        f.write(str(gene) + "\t" + str(gene_evaluation_median[gene]) + "\t" + str(gene_evaluation_full[gene]) + "\n")

    f.close()

if __name__ == "__main__":
    import os
    import sys

    print("Starting analysis...")

    #redirect stdout to nothing in order to speed up program
    old = sys.stdout
    null = open(os.devnull, 'w')
    #sys.stdout = null

    #run algo
    import timeit
    start = timeit.default_timer()
    run_analysis_on_dataset("BC", 10)
    end = timeit.default_timer()

    #reset output to terminal and print results!
    sys.stdout = old
    print("Took " + str(end - start) + "s")
