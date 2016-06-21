'''
Contains methods to fit genes to gaussian model as per the paper
'''

import timeit
import numpy as np

import matplotlib.pyplot as plt
from scipy.stats import norm
from db_setup import read_dumped_data
from sklearn.mixture import GMM

import gaussian_sampling

#for pickle formatting
from process_BC_data import sample
from process_BC_data import expression_profile
from process_BC_data import clinical_data

'''
x should data points of the intensity of the gene in question
y should be the corresponding class of each sample, supports only binary 1 and 0 for now

plot = boolean value whether to plot the resulting classifier's frequency curve
title = title of graph to plot

no test/train split since this is a clustering model

returns the trained model
'''
def fit_test_model(x):
    gauss_model = GMM(n_components=2, n_init=5, n_iter=1000)
    gauss_model.fit(x)

    return gauss_model

'''
Given a gauss mix model, prints parameters for each seperate peak
'''
def print_model_params(gauss_model):
    coeffs = gauss_model.weights_
    mus = [x[0] for x in gauss_model.means_]
    sigmas = [x[0] for x in gauss_model.covars_]

    string = ("Gaussian model: " + str(gauss_model)) + '\n'
    string += ("Coeff:\t" + str(coeffs)) + '\n'
    string += ("Mus:\t" + str(mus)) + '\n'
    string += ("Sigmas:\t" + str(sigmas)) + '\n'

    print(string)

'''
prints out information how the given model clusters and classified real data
Assuems 1D features
'''
def print_test_model(x, Y, gauss_model):
    Y_h = gauss_model.predict(x)

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

if __name__ == "__main__":
    genes = []
    for gene_sample in read_dumped_data("BC_gene_profiles.pkl"):
        genes.append(gene_sample)

    for gene in genes:
        if gene.name == "ERBB2":
            intensities = [[x] for x in gene.intensities]
            print(intensities)
            model = fit_test_model(intensities)
            #print_test_model(x, Y, model)
            print_model_params(model)

    #future list inc
    '''
    #first process data
    clinical_profiles = []
    for profile in read_dumped_data("BC_clinical_profiles.pkl"):
        clinical_profiles.append(profile)
    clinical_profiles.sort(key=lambda profile: profile.sample_num)

    samples = []
    for profile in read_dumped_data("BC_expression_profiles.pkl"):
       samples.append(profile)

    #extract training/testing data
    gene = "ERBB2"

    #ToDo: implement maping of num to profile and vice versa for that sexy linear time

    x = []
    Y = []
    for sample in samples:
        sample_num = sample.sample_num

        #extract the profile of interest
        profiles = sample.profiles
        gene_profile = ([p for p in profiles if p.gene == gene])[0]
        x.append([gene_profile.intensity])

        #check if patient has cancer
        for i in range(0, len(clinical_profiles)):
            if clinical_profiles[i].sample_num == sample_num:
                if clinical_profiles[i].posnodes == 'y' or \
                        clinical_profiles[i].event_meta == 1 or\
                        clinical_profiles[i].st_gallen == 0:
                    Y.append(1)
                else:
                    Y.append(0)
                break

    #train model
    model = fit_test_model(x)
    print_test_model(x, Y, model)
    print_model_params(model)
    '''



