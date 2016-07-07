'''
This file contains methods for importing and pre-processing the data from the Breast Cancer data set used in the
paper "Evaluating gene set enrichment analysis via a hybrid data model" by J. Hua et al. 2014.
'''

import os
import glob

from data_models import *

BC_DATA_DIR = os.getcwd() + "/Data/HybridSets/BC"
BC_EXPRESSION_DIR = BC_DATA_DIR + "/ExpressionProfiles"
BC_CLINICAL_DATA_FILE = BC_DATA_DIR + "/ClinicalData/ClinicalData.txt"
BC_GENE_SETS_FILE = BC_DATA_DIR + "/GeneSets/c2.all.v5.1.symbols.gmt"

class FileFormatError(Exception):
    """
    This does something
    """
    def __init__(self, file, error):
        self.file = file
        self.error = error

    def __str__(self):
        return "File " + str(self.file) + "\n" + str(self.error)

'''
Given a file containing expression profiles, reads and returns a list of correlating data_models
This assumes the gene expression profile format from the 295 sample study in the BC set. Note,
this format is encoded in utf-7!

file_name   =   the name of the file to read information from
'''
def readExpressionProfile(file_name):
    file = open(file_name, 'r', encoding="utf-7")

    #get rows of the file
    rows = file.read().split("\n")
    if(len(rows) < 2):
        raise FileFormatError(file, "Labels in top two rows missing")

    #extract sample id's from rows
    IDs = []
    first_row = rows[0].split("\t")
    first_row[:] = (val for val in first_row if val != '')
    for string in first_row:
        if "Sample" in string:
            IDs.append(int(string.split(' ')[1]))
        else:
            raise FileFormatError(file, "Sample IDs expected or incorrectly formatted")

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
        gene = fields[1]

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

'''
Reads all expression profile files in the given directory

dir     =       the directory containing the .txt files containing expression profiles to read from

Returns a list of expression profiles
'''
def getExpressionProfiles(dir):
    files = [filename for filename in glob.glob(os.path.join(dir, "*.txt"))]
    expression_profiles = []
    for file in files:
        expression_profiles += readExpressionProfile(file)

    return expression_profiles

'''
Given a file containing clinical profiles, reads and returns a list of correlating data_models
This assumes the clinical profile format from the 295 sample study in the BC set. Note,
this format is encoded in utf-7!

file_name   =   the name of the file to read information from
'''
def getClinicalData(file_name):
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

'''
Reads a file containing gene set data, one set per line, tab seperated where index 0 is the name
of the pathway, index 1 is the url for pathway info, and the rest the name of the genes

file_name   =   the name of the file to read information from
'''
def getGeneSetData(file_name):
    file = open(file_name, 'r', encoding="utf-7")

    gene_sets = []
    rows = file.read().split('\n')

    for row in rows:
        if(row == ''):
            break

        set_info = row.split("\t")
        set_name = set_info[0]
        set_url = set_info[1]

        genes = []
        for i in range(2, len(set_info)):
            genes.append(set_info[i])

        gene_sets.append(gene_set(set_name,set_url, genes))

    return gene_sets