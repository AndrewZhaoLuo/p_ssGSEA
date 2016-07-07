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

    print("Training model for every gene...")
    models = {}
    for gene in gene_profiles.keys():
        expressions = [[exp] for exp in gene_profiles[gene]]
        models[gene] = fit_test_model(expressions)

    genes = []
    for gene_sample in gene_profiles:
        genes.append(gene_sample)

    return models