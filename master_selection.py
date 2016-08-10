'''
This file contains a series of scripts for evaluating gene models and determining their usability as
a "master gene."
'''
import math

from scipy.stats import norm

def calculate_prior(model):
    """
    Given a gaussian mixture model fitted to the expression of a single gene, returns the prior calculation of the gene

    :requires: the model is composed of two gaussian components!

    :param model: the mixture model of a given gene
    :type model: sklearn.mixture.GMM

    :returns: float representing the calculated prior of this model
    """
    coeffs = model.weights_

    return min(coeffs[0], 1 - coeffs[0])

def solveQuadratic(a, b, c):
    """
    Given coefficients a, b, and c. Solves the quadratic in the form :math:`ax^2 + bx + c = 0`.

    :requires: the determinant of the quadratic is non-zero

    :param a: the coefficient of the x^2 term
    :type a: float

    :param b: the coefficient of the x term
    :type b: float

    :param c: the coefficient of 1 term
    :type c: float

    :returns: a tuple with two values containing two roots of the equation. If there is only one distinct root,
              it will appear twice
    """
    det = b ** 2 - 4 * a * c
    r1 = (-b + det ** 0.5) / (2 * a)
    r2 = (-b - det ** 0.5) / (2 * a)

    return (r1, r2)

'''
Given the sigma and mu of two distributions, returns the x value where the two distributions are equal
'''
def findIntersection(mu1, sigma1, mu2, sigma2):
    """
    Finds and returns the intersection points of two standard deviations

    :requires: the parameters given describe two different standard deviations

    :param mu1: the mean of the first distribution
    :type mu1: float

    :param sigma1: the standard deviation of the first distribution
    :type sigma1: float

    :param mu2: the mean of the second distribution
    :type mu2: float

    :param sigma2: the standard deviation of the second distribution
    :type sigma2: float

    :returns: the intersection points of the two distributions as a tuple pair. If there is only one intersection,
              the point appears twice in the tuple.
    """

    #find coefficient of quadratic
    a = sigma2 ** 2 - sigma1 ** 2
    b = 2 * sigma1 ** 2 * mu2 - 2 * sigma2 ** 2 * mu1
    c = mu1 ** 2 * sigma2 ** 2 - mu2 ** 2 * sigma1 ** 2 - 2 * sigma1 ** 2 * sigma2 ** 2 * math.log(sigma2 / sigma1)

    #find roots of the quadratic
    return solveQuadratic(a, b, c)

'''
Given two distributions, returns the x coordinate of the decision boundary
based on where both have equal z scores
'''
def findDecisionBoundary(mu1, sigma1, mu2, sigma2):
    """
    Given the parameters describing two gaussian models representing two classes, returns the decision boundary
    for classification. THis assumes classification is done through a *z-score*.

    :requires: the parameters given describe two different standard deviations

    :param mu1: the mean of the first distribution
    :type mu1: float

    :param sigma1: the standard deviation of the first distribution
    :type sigma1: float

    :param mu2: the mean of the second distribution
    :type mu2: float

    :param sigma2: the standard deviation of the second distribution
    :type sigma2: float

    :returns: the decision boundary as a float

    """
    return (mu1 * sigma2 - mu2 * sigma1) / (sigma2 - sigma1)

'''
Returns the bayes error of the given mixture model
'''
def calculate_bayes_error(model):
    """
    Returns the bayes error of the given mixture model. This is calculated by integrating the false class density
    around the decision boundary of the model.

    :requires: the model is composed of two gaussian components!

    :param model: the mixture model of a given gene
    :type model: sklearn.mixture.GMM

    :returns: the bayes error of the classifier as a float
    """
    #first we find the intersection point
    coeffs = model.weights_
    mus = [x[0] for x in model.means_]
    sigmas = [x[0] ** 0.5 for x in model.covars_]
    r1, r2 = findIntersection(mus[0], sigmas[0], mus[1], sigmas[1])

    root = 0
    if r1 < max(mus[0], mus[1]) and r1 > min(mus[0], mus[1]):
        root = r1
    else:
        root = r2

    #now that we have the intersectionm we need the CDF/survival function of both plots
    err = 0
    if(root < mus[0]):
        err += norm.sf(root, loc=mus[1], scale=sigmas[1]) * coeffs[1]
        err += norm.cdf(root, loc=mus[0], scale=sigmas[0]) * coeffs[0]
    else:
        err += norm.sf(root, loc=mus[0], scale=sigmas[0]) * coeffs[0]
        err += norm.cdf(root, loc=mus[1], scale=sigmas[1]) * coeffs[1]

    return err #/ (norm.sf(-10000, loc=mus[0], scale=sigmas[0]) + norm.sf(-10000, loc=mus[1], scale=sigmas[1]) - err)

def calculate_populairity(gene, gene_sets):
    '''
    Returns the popularity of the given gene with respect to the given gene_sets. This is the number of times the
    gene appears in the given gene_sets

    :param gene: the gene to check the popularity of
    :type gene: str

    :param gene_sets: a map of gene_set names to gene_set objects
    :type gene_sets: dict

    :returns: the number of times the gene appears in the given gene sets
    '''

    return sum([1 for key in gene_sets.keys() if gene in gene_sets[key].genes])

def calculate_fold_change(model):
    """
    Returns the fold change of the given mixture model

    :requires: the model is composed of two gaussian components!

    :param model: the mixture model of a given gene
    :type model: sklearn.mixture.GMM

    :returns: the fold change of the given gene model as a float
    """
    mus = model.means_
    return abs(mus[0] - mus[1])

def calculate_shape_balance(model):
    """
    Returns the shape imbalance of the given model

    :requires: the model is composed of two gaussian components!

    :param model: the mixture model of a given gene
    :type model: sklearn.mixture.GMM

    :returns: The shape imbalance of the gene model
    """
    sigmas = model.covars_

    return max(sigmas[1] / sigmas[0], sigmas[0] / sigmas[1])