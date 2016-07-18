'''
This file contains methods for plotting various graphs given various data points
'''

import os
import pickle

a = pickle.load(open(os.getcwd() + '/Data/AppCache/BC/CachedEnrichmentPValueSplit.pkl', 'rb'))

for gene in range[0]