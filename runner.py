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
from cache_codec import load_bayes_scores
from cache_codec import load_null_scores

from analyze_enrichment import rank_by_t_test_keyed
from analyze_enrichment import evaluate_rankings_keyed

import pickle
import random

DEFAULT_NUM_PROCESSES = 7

def run_analysis_on_dataset(NUM_PROCESSES, data_set, n, pheno_sample, gene_options='all', test='ssGSEA'):
    '''
    Given data_set, and n which is the number of replicates for the experiment:
    for every gene in the data_set for which there exists a filtered gene set, creates imaginary phenotypes for every
    sample using a gaussian mixture model. Imaginary phenotypes are then used to fit gaussian models of enrichment values
    for each gene set. A gene set is considered "intersting" if the generated models can seperate. This is done
    by a t-test of the class sepearated enrichment data. Gene ranking is done by looking at the sets in which the
    master gene belongs. The master gene should show up at the very top of the sets.

    ToDo: finish doc

    ToDo: get rid of hackey things done in this evil abomination of code. My mentors are crying in their sleep
    '''

    print("Starting analysis...")
    master_genes = load_gene_popularity(data_set).keys()

    enrichment_scores = None
    if test == 'ssGSEA':
        enrichment_scores = load_ssGSEA_scores(data_set)
    elif test == 'bayes_low':
        enrichment_scores = load_bayes_scores(data_set, 'low')
    elif test == 'bayes_mid':
        enrichment_scores = load_bayes_scores(data_set, 'mid')
    elif test == 'bayes_high':
        enrichment_scores = load_bayes_scores(data_set, 'high')
    elif test == "null":
        enrichment_scores = load_null_scores(data_set)
    else:
        print("ERROR, invalid enrichment test!")
        return

    gene_sets = load_filtered_gene_sets(data_set)
    best_models = load_best_models("BC", 10, 10)

    good_genes = []
    for master_gene in master_genes:
        if ((master_gene in best_models.keys() and gene_options == 'all') or master_gene in gene_options) and master_gene != '':
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

    #dataman!!
    from cache_codec import DATA_DIR
    pickle.dump(enrichment_map, open(DATA_DIR + "/" + data_set + "/" + str(test) + "CachedEnrichmentPValueSplit.pkl", 'wb'))

    gene_evaluation_median = {}
    gene_evaluation_full = {}

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
    f = open(test + 'results_rankings.txt','w')
    f.write("Gene\tMedian Rank\tFull Ranks\n")
    for gene in sorted_genes:
        f.write(str(gene) + "\t" + str(gene_evaluation_median[gene]) + "\t" + str(gene_evaluation_full[gene]) + "\n")
    f.close()

def main(argv):
    if len(argv) < 4:
        print("Usage: python runner.py [NUM_PROCESSES] [REPLICATES] [SAMPLES_PER_REPLICATE] [TEST: ssGSEA | bayes_high | bayes_low | bayes_mid] <...GENES, empty for all>")
        return

    NUM_PROCESSES = int(argv[0])
    REPLICATES = int(argv[1])
    SAMPLES_PER_REPLICATE = int(argv[2])
    TEST = str(argv[3])

    genes = []
    if len(argv) < 5:
        genes = 'all'
    else:
        genes = [gene for gene in argv[4:]]

    print("Starting analysis...")

    #run algo
    import timeit
    start = timeit.default_timer()
    run_analysis_on_dataset(NUM_PROCESSES, "BC", REPLICATES, SAMPLES_PER_REPLICATE, gene_options=genes, test=TEST)
    end = timeit.default_timer()

    #reset output to terminal and print results!
    print("Took " + str(end - start) + "s")

if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
