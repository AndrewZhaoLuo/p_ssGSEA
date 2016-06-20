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

print("Processing data...")
start = timeit.default_timer()

#first process data
clinical_profiles = []
for profile in read_dumped_data("BC_clinical_profiles.pkl"):
    clinical_profiles.append(profile)

clinical_profiles.sort(key=lambda profile: profile.sample_num)

samples = []
for profile in read_dumped_data("BC_expression_profiles.pkl"):
   samples.append(profile)

end = timeit.default_timer()
print("\tDone! Took" + str(end - start) + "s")

print("Grabbing genes...")
start = timeit.default_timer()


#extract gene of interest from people
#splitting into two groups: those with pos nodes and those without (cancer)
gene = "ERBB2"

'''
ToDo: implement maping of num to profile and vice versa for that sexy linear time
'''
pos_expression_profiles = []
neg_expression_profiles = []
for sample in samples:
    sample_num = sample.sample_num

    #extract the profile of interest
    profiles = sample.profiles
    gene_profile = ([p for p in profiles if p.gene == gene])[0]

    #check if patient has cancer
    for i in range(0, len(clinical_profiles)):
        if clinical_profiles[i].sample_num == sample_num:
            if clinical_profiles[i].st_gallen == 0:
                pos_expression_profiles.append(gene_profile)
            else:
                neg_expression_profiles.append(gene_profile)
            break

end = timeit.default_timer()
print("\tDone! Took" + str(end - start) + "s")

print("Fitting Models...")
start = timeit.default_timer()

pos_expression_intensity = [[p.intensity] for p in pos_expression_profiles]
neg_expression_intensity = [[p.intensity] for p in neg_expression_profiles]

gauss_both = GMM(n_components=2, n_init=5, n_iter=1000)
gauss_both.fit(pos_expression_intensity + neg_expression_intensity)
coeffs = gauss_both.weights_
mus = [x[0] for x in gauss_both.means_]
sigmas = [x[0] for x in gauss_both.covars_]

print("Coeff:\t" + str(coeffs))
print("Mus:\t" + str(mus))
print("Sigmas:\t" + str(sigmas))

samples = gauss_both.sample(100000)
gaussian_sampling.plot_multidist("ERBB2", samples)


end = timeit.default_timer()
print("\tDone! Took" + str(end - start) + "s")


