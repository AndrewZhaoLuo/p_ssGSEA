'''
This file contains methods for importing and pre-processing the data from the Breast Cancer data set used in the
paper "Evaluating gene set enrichment analysis via a hybrid data model" by J. Hua et al. 2014.

ASSUMES ALL FILES ARE ENCODED AS UTF-7! (Which they are originally!)
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

def getClinicalData(file_name):
    """
    Reads clinical information about each patient who gave a gene sample

    :param file_name: the directory and file containing the patient's clinical information
    :type file_name: str

    :return: a list of clinical_data
    """
    file = open(file_name, 'r', encoding="utf-7")

    clinical_datas = []
    rows = file.read().split('\n')

    #ignore first line of labels
    #has extra blank line so don't take last row
    for i in range(1, len(rows) - 1):
        row = rows[i].split('\t')

        sample_num      = row[0]

        first_series    = row[1]
        posnodes        = row[2]
        event_meta      = row[3]
        event_death     = row[4]
        time_survival   = row[5]
        time_recur      = row[6]
        time_meta       = row[7]
        esr1            = row[8]
        nih             = row[9]
        st_gallen       = row[10]
        conserv         = row[11]
        c1_from_data    = row[12]
        c1_cross_valid  = row[13]
        c1_used         = row[14]

        other_fields = {}
        other_fields["first_series"] = first_series
        other_fields["posnodes"] = posnodes
        other_fields["event_meta"] = event_meta
        other_fields["event_meta"] = event_meta
        other_fields["event_death"] = event_death
        other_fields["time_survival"] = time_survival
        other_fields["time_recur"] = time_recur
        other_fields["time_meta"] = time_meta
        other_fields["esr1"] = esr1
        other_fields["nih"] = nih
        other_fields["st_gallen"] = st_gallen
        other_fields["conserv"] = conserv
        other_fields["c1_from_data"] = c1_from_data
        other_fields["c1_cross_valid"] = c1_cross_valid
        other_fields["c1_used"] = c1_used

        data = clinical_data(sample_num, other_fields)
        clinical_datas.append(data)

    return clinical_datas