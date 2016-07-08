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
def simulate_data(models, master_gene_name, sample_profiles, n):
    """
    Given a list of gaussian mixture models for genes, a master gene, and a list of sample profiles,
    returns n simulated phenotypes for each sample profile. This is done by selection of sample profiles
    without replacement n times. Note this and the actual simulation is a stochastic process which looks
    at the combined likehood for each class according to the mixture model.

    :requires: n is positive and less than or equal to len(sample_profiles)

    :param models: a dictionary of strings as gene names to the gaussian mixture model (GMM) for each gene.
    :type models: dict

    :param master_gene_name: the gene to be used as the maste gene for phenotype simulation
    :type master_gene_name: str

    :param sample_profiles: a list of samples (data_models.sample)
    :type sample_profiles: list

    :param n: the number of phenotypes to simulate
    :type n: int

    :return: a dict mapping sample id's to phenotype classes. This phenotype does not necessarily correspond
             to any real phenotype but rather the components of the mixture model of the master gene. class 0 represents
             the distribution with a smaller mean while class 1 represents the one with the larger mean.
    """

    master_gene_model = models[master_gene_name]

    #for selecting random front n values
    random.shuffle(sample_profiles)

    classifications = {}
    for i in range(0, n):
        cur_profile = sample_profiles[i]
        id = cur_profile.id
        cur_intensity = cur_profile.profiles[master_gene_name].intensity

        #classification of phenotype for this profile
        proportions = master_gene_model.weights_
        mus = master_gene_model.means_
        sigmas = [x[0] ** 0.5 for x in master_gene_model.covars_]

        #setting 0/1 to the proper low/high peaks. 0 = smaller mu
        if(mus[0] < mus[1]):
            mu0 = sum(mus[0])
            sigma0 = sigmas[0]
            proportions0 = proportions[0]
            mu1 = sum(mus[1])
            sigma1 = sigmas[1]
            proportions1 = proportions[1]
        else:
            mu0 = sum(mus[1])
            sigma0 = sigmas[1]
            proportions0 = proportions[1]
            mu1 = sum(mus[0])
            sigma1 = sigmas[0]
            proportions1 = proportions[0]

        #figuring out probability for each sample being in each peak
        density0 = norm.pdf(cur_intensity, loc=mu0, scale=sigma0)
        density1 = norm.pdf(cur_intensity, loc=mu1, scale=sigma1)

        probability1 = proportions1 * density1 / (proportions1 * density1 + proportions0 * density0)

        if random.random() <= probability1:
            classifications[id] = 1
        else:
            classifications[id] = 0

    return classifications
if __name__ == "__main__":
    data_sets = pickle.load(open("BC_simulated_phenotypes2.pkl", 'rb'))

    #maps id -> number of times phenotype 0/1 appear
    count0 = {}
    count1 = {}

    #initializae map
    for key in data_sets[0].keys():
        count0[key] = 0
        count1[key] = 0

    for i in range(0, len(data_sets)):
        set = data_sets[i]
        for key in set.keys():
            if set[key] == 0:
                count0[key] += 1
            else:
                count1[key] += 1

    #find a sort by class imbalance
    diffRatio = {}
    for key in data_sets[0].keys():
        diffRatio[key] = abs(count0[key] - count1[key])

    sorted_keys = sorted(diffRatio, key=diffRatio.get)
    print("key\t\tclass_0\tclass_1")
    for key in sorted_keys:
        print(str(key) + "\t\t" + str(count0[key]) + "\t\t" + str(count1[key]))


