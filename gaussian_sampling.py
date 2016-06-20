'''
This file contains code snippets demonstrating one way
to create and sample points from a gaussian distribution,
plot those points, and see how definition of the curve changes as
the number of samples of the curve is taken.
'''

from random import gauss
import matplotlib.pyplot as plotter
import seaborn as sea

SIGMA = [100, 1000, 10000]
SAMPLE_NUMBERS = [10, 100, 1000]

'''
Returns a set of points sampled from a distribution

num_samples     =       number of points to sample, must be postive
mu              =       mean of the dist.
sigma           =       standard deviation of the dist, must be positive

returns an array of size num_samples containing points sampled from the dist.
'''
def sample_dist(num_sample, mu, sigma, coeff):
    assert num_sample > 0
    assert sigma > 0

    return [gauss(mu, sigma) * coeff for i in range(0, num_sample)]

'''
Returns a set of points sampled from multiple gaussian dist!

as sample_dist except takes all lists/arrays. the arrays should
all be the same length where parameters sharing the same index
will belong to one distribution

returns all the sampling of all distributions as one list
'''
def sample_multi_dist(num_samples, mus, sigmas, coeffs):
    assert len(num_samples) == len(mus)
    assert len(mus) == len(sigmas)

    samples = []
    for i in range(0, len(num_samples)):
        samples += sample_dist(num_samples[i], mus[i], sigmas[i], coeffs[i])

    return samples

'''
plots graphs of distributions as sampling and sigma changes and exports it as a png
'''
def plot_sigma_samplesize(sample_numbers, sigmas, size, title):
    # partition subplots
    f, axes = plotter.subplots(len(sample_numbers), len(sigmas), figsize=size, sharex= False, sharey= False)
    #sea.despine(left=True)

    for a in range(0, len(sample_numbers)):
        for b in range(0, len(sigmas)):
            sample_num = sample_numbers[a]
            sigma = sigmas[b]

            values = sample_dist(sample_num, mu=0, sigma=sigma)
            sea.distplot(values, ax=axes[b,a], axlabel=("#SAMPLES=" + str(sample_num) + "    " + "SIGMA=" + str(sigma)))

    plotter.savefig(title)
    plotter.close()

'''
plots graph of sampling from multiple distributions and exports it as a png

see sample_multi_dist for instructions on parameters
'''
def plot_multidist(num_samples, mus, sigmas, coeffs, title, values):
    values = sample_multi_dist(num_samples, mus, sigmas, coeffs)
    sea.distplot(values)

    plotter.title("SAMPLE=" + str(num_samples) + " MU=" + str(mus) + " SIGMA=" + str(sigmas))
    plotter.savefig(title)
    plotter.close()

def plot_multidist(title, values):
    sea.distplot(values)

    plotter.title(title)
    plotter.savefig(title)
    plotter.close()


if __name__ == "__main__":
    #plot_sigma_samplesize(SAMPLE_NUMBERS, SIGMA, (24,24), "Single_Distribution_Test")
    plot_multidist([10000, 10000], [0, 0], [100, 500], [1,1], "Double_Distribution_Test")
