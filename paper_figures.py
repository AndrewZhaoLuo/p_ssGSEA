'''
Figures for the poster
'''

import matplotlib

from cache_codec import *
from matplotlib import pyplot as plotter

import seaborn

if __name__ == "__main__":
	#Figure showing gaussian mixture model for ERBB2
	models = load_gene_models("BC")
	samples = load_sample_profiles("BC").values()
	expression_levels = [sample.profiles["ERBB2"].intensity for sample in samples]

	print(expression_levels)

	seaborn.rugplot(expression_levels, color="black")
	seaborn.distplot(models["ERBB2"].sample(100000), hist=False)

	exp_patch = matplotlib.patches.Patch(color='black', label="Expression Data")
	gmm_patch = matplotlib.patches.Patch(color='blue', label="GMM model")
	plotter.legend(handles=[exp_patch, gmm_patch])
	plotter.legend(loc="upper right")

	plotter.title("ERBB2")
	plotter.xlabel("Expression Level")
	plotter.ylabel("Density")

	plotter.savefig("Yay!.png")
