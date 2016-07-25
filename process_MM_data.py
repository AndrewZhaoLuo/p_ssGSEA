'''
Scripts for dealing with importing MM dataset
'''

import os
import glob

from data_models import *

BC_DATA_DIR = os.getcwd() + "/Data/HybridSets/BC"
BC_EXPRESSION_DIR = BC_DATA_DIR + "/ExpressionProfiles"
BC_CLINICAL_DATA_FILE = BC_DATA_DIR + "/ClinicalData/ClinicalData.txt"


class BC_FileFormatError(Exception):
    """
    This class represents an error in formatting from a file from the BC data set.
    """
    def __init__(self, file, error):
        self.file = file
        self.error = error

    def __str__(self):
        return "File " + str(self.file) + "\n" + str(self.error)

def readExpressionProfile(file_name):
    """
    Reads a single file from the BC set containing gene expression information.

    :param file_name: the path and file containing expression information to read
    :type file_name: str

    :return: a list of expression_profiles
    """
    file = open(file_name, 'r', encoding="utf-7")

    #get rows of the file
    rows = file.read().split("\n")
    if(len(rows) < 2):
        raise BC_FileFormatError(file, "Labels in top two rows missing")

    #extract sample id's from rows
    IDs = []
    first_row = rows[0].split("\t")
    first_row[:] = (val for val in first_row if val != '')
    for string in first_row:
        if "Sample" in string:
            IDs.append(int(string.split(' ')[1]))
        else:
            raise BC_FileFormatError(file, "Sample IDs expected or incorrectly formatted")

    #extract expression data into here
    expression_profiles = []

    #start from 2 which is the first row that contains gene expression data
    for i in range(2, len(rows)):
        row = rows[i]

        #some of the data has a random newline
        if str(row) == '':
            break

        #each row contains multiple samples for each gene45
        fields = row.split("\t")

        substance = fields[0]

        #sanitize / to prevent issues with directories being created
        gene = fields[1].replace('/', '-')

        #start from two, which is the index of the first column with actual information
        #5 fields, so offset by 5 with each id
        for f in range(0, len(IDs)):
            sample_num  = IDs[f]
            log_ratio   = fields[2 + 5 * f + 0]
            log_error   = fields[2 + 5 * f + 1]
            p_value     = fields[2 + 5 * f + 2]
            intensity   = fields[2 + 5 * f + 3]
            flag        = fields[2 + 5 * f + 4]

            other_fields = {}
            other_fields["log_ratio"] = log_ratio
            other_fields["log_error"] = log_error
            other_fields["p_value"] = p_value
            other_fields["flag"] = flag
            other_fields["substance"] = substance

            profile = expression_profile(sample_num, gene, intensity, other_fields)
            expression_profiles.append(profile)

    file.close()
    return expression_profiles


def readAllExpressionProfiles(dir):
    """
    Reads all files containing expression data in the given directory

    :param dir: the directory from where to read all text files containing expression data.
    :type dir: str

    :return: a list of expression_profiles formed from all files read in the directory
    """
    files = [filename for filename in glob.glob(os.path.join(dir, "*.txt"))]
    expression_profiles = []
    for file in files:
        expression_profiles += readExpressionProfile(file)

    return expression_profiles
