'''
This module contains methods for simulating data given a hybrid data model
'''
import pickle
import random
from scipy.stats import norm

'''
Given a master_gene's model and a list of sample profiles, returns a
map of sample_nums -> phenotype classifications for this set.

0 = peak with peak lower on x axis
1 = peak with peak higher on x axis

n samples specifically
'''
def simulate_data(master_gene_model, master_gene_name, sample_profiles, n):
    #for selecting random front n values
    random.shuffle(sample_profiles)

    classifications = {}
    for i in range(0, n):
        cur_profile = sample_profiles[i]
        sample_num = cur_profile.sample_num
        cur_intensity = cur_profile.profiles[master_gene_name].intensity

        #classification of phenotype for this profile
        proportions = master_gene_model.weights_
        mus = master_gene_model.means_
        sigmas = [x[0] ** 0.5 for x in master_gene_model.covars_]

        #1 = peak with lesser mu, 2 = peak with greater mu
        mu1 = 0
        sigma1 = 0
        proportions1 = 0
        mu2 = 0
        sigma2 = 0
        proportions2 = 0

        #setting 1/2 to the proper low/high peaks
        if(mus[0] < mus[1]):
            mu1 = sum(mus[0])
            sigma1 = sigmas[0]
            proportions1 = proportions[0]
            mu2 = sum(mus[1])
            sigma2 = sigmas[1]
            proportions2 = proportions[1]
        else:
            mu1 = sum(mus[1])
            sigma1 = sigmas[1]
            proportions1 = proportions[1]
            mu2 = sum(mus[0])
            sigma2 = sigmas[0]
            proportions2 = proportions[0]

        #figuring out probability for each sample being in each peak
        density1 = norm.pdf(cur_intensity, loc=mu1, scale=sigma1)
        density2 = norm.pdf(cur_intensity, loc=mu2, scale=sigma2)

        probability1 = proportions1 * density1 / (proportions1 * density1 + proportions2 * density2)

        if random.random() <= probability1:
            classifications[sample_num] = 0
        else:
            classifications[sample_num] = 1

    return classifications

if __name__ == "__main__":
    master_genes = pickle.load(open("BC_master_genes.pkl", 'rb'))
    models = pickle.load(open("BC_trained_models.pkl", 'rb'))
    profiles = pickle.load(open("BC_expression_profiles.pkl", 'rb'))

    ERBB2_model = models["ERBB2"]
    g = simulate_data(ERBB2_model, "ERBB2", profiles, 100)
    print(g)