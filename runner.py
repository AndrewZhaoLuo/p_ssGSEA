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
from cache_codec import load_best_models

from analyze_enrichment import rank_by_t_test_keyed
from analyze_enrichment import evaluate_rankings_keyed

import numpy.random as nrandom
import random

NUM_PROCESSES = 8

def run_analysis_on_dataset(data_set, n, pheno_sample, gene_options='all'):
    '''
    Given data_set, and n which is the number of replicates for the experiment:
    for every gene in the data_set for which there exists a filtered gene set, creates imaginary phenotypes for every
    sample using a gaussian mixture model. Imaginary phenotypes are then used to fit gaussian models of enrichment values
    for each gene set. A gene set is considered "intersting" if the generated models can seperate. This is done
    by a t-test of the class sepearated enrichment data. Gene ranking is done by looking at the sets in which the
    master gene belongs. The master gene should show up at the very top of the sets.

    ToDo: finish doc
    '''

    print("Starting analysis...")
    master_genes = load_gene_popularity(data_set).keys()
    enrichment_scores = load_ssGSEA_scores(data_set)
    gene_sets = load_filtered_gene_sets(data_set)
    best_models = load_best_models("BC", 10, 10) 

    gene_evaluation_median = {}
    gene_evaluation_full = {}

    good_genes = []
    for master_gene in master_genes:
        #if gene_options == 'all' or master_gene in gene_options:
        if master_gene in best_models.keys()
            good_genes.append(master_gene)
    print("Loaded data!")

    #ToDo: find way to do this without spawning proccesses
    pool = Pool(processes=NUM_PROCESSES)

    #load phenotypes, maps master gene to phenotype
    pheno_params = [(data_set, n, good_gene) for good_gene in good_genes]
    loaded_phenos = pool.starmap(load_sim_phenotype_keyed, pheno_params)

    #process into more usable form of master gene -> phenotypes
    pheno_map = {}
    for map in loaded_phenos:
        for key in map.keys():
            sampled_pheno_map = []
            for phenotype in map[key]:
                sampled_keys = random.sample(phenotype.keys(), pheno_sample)
                sampled_map = {id: phenotype[id] for id in sampled_keys}
                sampled_pheno_map.append(sampled_map)

            pheno_map[key] = sampled_pheno_map
    print("Calculated phenotypes")

    #work out the enrichment_score ranks
    enrichment_list = pool.starmap(rank_by_t_test_keyed, [(enrichment_scores, pheno_map[gene], gene) for gene in good_genes])
    enrichment_map = {}
    for enrichment in enrichment_list:
        for key in enrichment.keys():
            enrichment_map[key] = enrichment[key]
    print("Calculated enrichment")

    #finally translate that to gene rankings
    results_map = pool.starmap(evaluate_rankings_keyed, [(enrichment_map[gene], gene_sets, gene) for gene in good_genes])
    for result in results_map:
        for master_gene in result.keys():
            evaluation = result[master_gene]
            if len(evaluation) > 0:
                print(master_gene + "\t" + str(evaluation))
                gene_evaluation_median[master_gene] = numpy.median(evaluation)
                gene_evaluation_full[master_gene] = evaluation
    print("Calculated ranking")

    pool.close()
    pool.join()

    #print out final rankings
    sorted_genes = sorted(gene_evaluation_median, key=gene_evaluation_median.get)
    f = open('results.txt','w')
    f.write("Gene\tMedian Rank\tFull Ranks\n")
    for gene in sorted_genes:
        f.write(str(gene) + "\t" + str(gene_evaluation_median[gene]) + "\t" + str(gene_evaluation_full[gene]) + "\n")

    f.close()

import os
import sys

def main(argv):
    if len(argv) < 3:
        print("Usage: python runner.py [NUM_PROCESSES] [REPLICATES] [SAMPLES_PER_REPLICATE] <...GENES>")
        return

    NUM_PROCESSES = int(argv[0])
    REPLICATES = int(argv[1])
    SAMPLES_PER_REPLICATE = int(argv[2])

    genes = []
    if len(argv) < 4:
        genes = 'all'
    else:
        genes = [gene for gene in argv[3:]]

    print("Starting analysis...")

    #run algo
    import timeit
    start = timeit.default_timer()
    run_analysis_on_dataset("BC", REPLICATES, SAMPLES_PER_REPLICATE, gene_options=genes)
    end = timeit.default_timer()

    #reset output to terminal and print results!
    print("Took " + str(end - start) + "s")

if __name__ == "__main__":
    main(sys.argv[1:])
