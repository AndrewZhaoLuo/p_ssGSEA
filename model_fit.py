'''
Contains methods to fit gene expression to mixed gaussian model as per the paper
'''

import timeit
import numpy as np

import matplotlib.pyplot as plt
from scipy.stats import norm
from db_setup import read_dumped_data
from sklearn.mixture import GMM

import gaussian_sampling
import pickle

#for pickle formatting
from process_BC_data import sample
from process_BC_data import expression_profile
from process_BC_data import clinical_data
from process_BC_data import gene_set

'''
Given an array of data points returns a 2 component gaussian mixture model which best estimates
the distribution of values. Uses E-M algorithm to accomplish this task

x       =       the data to cluster around

returns the trained model
'''
def fit_test_model(x):
    model = GMM(n_components=2, n_init=5, n_iter=10000, covariance_type='diag')
    model.fit(x)

    return model

'''
Given a gauss. mix model, prints parameters for each gaussian component
'''
def print_model_params(gauss_model):
    coeffs = gauss_model.weights_
    mus = [x[0] for x in gauss_model.means_]
    sigmas = [x[0] ** 0.5 for x in gauss_model.covars_]

    string = ("Gaussian model: " + str(gauss_model)) + '\n'
    string += ("Coeff:\t" + str(coeffs)) + '\n'
    string += ("Mus:\t" + str(mus)) + '\n'
    string += ("Sigmas:\t" + str(sigmas) + '\n')

    print(string)

'''
Evalautes the predictive power of a model and prints out confusion matrix of results

x       =       an array of expression profiles for a gene
Y       =       an array of phenotype classifications. should be done so the nth index of Y corresponds to
                    the nth index of array x
model   =       the model to evluate
'''
def evaluate_model(x, Y, model):
    Y_h = model.predict(x)

    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for i in range(0, len(Y)):
        if Y[i] == 1 and Y_h[i] == 1:
            tp += 1
        elif Y[i] == 1 and Y_h[i] == 0:
            fn += 1
        elif Y[i] == 0 and Y_h[i] == 0:
            tn += 1
        elif Y[i] == 0 and Y_h[i] == 1:
            fp += 1

    tab = "\t\t\t\t"
    print(tab + tab + "Actual")
    print(tab + tab + "Positive" + tab + "Negative")
    print("Predicted Positive" + tab + str(tp) + tab + str(fp))
    print("Predicted Negative" + tab + str(fn) + tab + str(tn))