'''
Contains methods to fit gene expression to mixed gaussian model using the Expectation Maximization algorithm
'''

from sklearn.mixture import GMM

def fit_test_model(x):
    """
    Given a series of samples, representing gene expression, returns a 2-component gaussian mixture model fitted
    to the samples. Expectation-Maximization is used to fit the data

    :param x: a list of floats representing gene expression profile samples for one gene. each entry in the list
              should be vectorized! ie. [[1.0], [3.2], [4.3]]
    :type x: list

    :returns: a sklearn.mixture.GMM, which represents a two component gaussian mixture model
    """
    model = GMM(n_components=2, n_init=5, n_iter=10000, covariance_type='diag')
    model.fit(x)

    return model

def print_model_params(gauss_model):
    """
    Prints out the coeff, mean, and standard deviation of the given mixture model. Also prints out information
    related to how the model was built

    :param model: the mixture model of a given gene
    :type model: sklearn.mixture.GMM
    """
    coeffs = gauss_model.weights_
    mus = [x[0] for x in gauss_model.means_]
    sigmas = [x[0] ** 0.5 for x in gauss_model.covars_]

    string = ("Gaussian model: " + str(gauss_model)) + '\n'
    string += ("Coeff:\t" + str(coeffs)) + '\n'
    string += ("Mus:\t" + str(mus)) + '\n'
    string += ("Sigmas:\t" + str(sigmas) + '\n')

    print(string)

def get_trained_models(gene_profiles):
    """
    Given a mapping of genes to expression values, returns a gaussian model for each gene

    :param gene_profiles: a dictionary of strings representing gene names to lists of floats representing
                          expression values
    :type gene_profiles: dict

    :return: a dict mapping gene names to the trained model of the gene
    """

    print("\tTraining model for every gene...")
    i = 0
    models = {}
    for gene in gene_profiles.keys():
        expressions = [[exp] for exp in gene_profiles[gene]]
        models[gene] = fit_test_model(expressions)

        i += 1
        if i % 100 == 0:
            print("Trained " + str(i) + " out of " + str(len(gene_profiles.keys())))

    genes = []
    for gene_sample in gene_profiles:
        genes.append(gene_sample)

    return models

if __name__ == "__main__":

    '''
    import master_selection

    file = open("./6-12-16Analysis/genes", 'r')
    genes = file.read().split('\n')

    out = open("pops", 'w')
    for gene in genes:
        pop = master_selection.calculate_populairity(gene, cache_codec.load_filtered_gene_sets("BC"))
        out.write(gene + "\t" + str(pop) + "\n")

    for set in cache_codec.load_filtered_gene_sets("BC"):
        if "ERBB2" in cache_codec.load_filtered_gene_sets("BC")[set].genes:
            print("\t" + set)
        if "SMID_BREAST_CANCER_ERBB2_DN" == set:
            print(cache_codec.load_filtered_gene_sets("BC")[set].genes)
    '''

    import cache_codec
    #code given master gene, prints out graph of
    #generated class0/1 phenotype profiles
    MASTER_GENE = "CLIC3"

    gene_sets = cache_codec.load_filtered_gene_sets("BC")
    enrichment_scores = cache_codec.load_ssGSEA_scores("BC")

    phenotypes = cache_codec.load_sim_phenotypes("BC", 10, MASTER_GENE)

    model = cache_codec.load_gene_models("BC")[MASTER_GENE]
    good_sets = [sets for sets in enrichment_scores.keys() if MASTER_GENE in gene_sets[sets].genes]

    import Extra_Modules.gaussian_sampling
    import analyze_enrichment
    for set in good_sets:
        class0, class1 = analyze_enrichment.analyze_phenotype_score_dist(enrichment_scores, phenotypes[0], set)
        print(class0)
        print(class1)
        Extra_Modules.gaussian_sampling.plot_multidist_from_values(MASTER_GENE + "_" + set,[class0, class1])

