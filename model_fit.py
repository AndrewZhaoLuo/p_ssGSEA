'''
Contains methods to fit genes to gaussian model as per the paper
'''

import timeit

from process_BC_data import sample
from process_BC_data import expression_profile
from process_BC_data import clinical_data

from db_setup import read_dumped_data
import pickle

start = timeit.default_timer()

expression_profiles = []
for profile in read_dumped_data("BC_expression_profiles.pkl"):
   expression_profiles.append(profile)

clinical_profiles = []
for profile in read_dumped_data("BC_clinical_profiles.pkl"):
    clinical_profiles.append(profile)

end = timeit.default_timer()

print((end - start))
