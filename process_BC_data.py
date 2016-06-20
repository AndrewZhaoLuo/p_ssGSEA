'''
This file contains methods for importing and pre-processing the sample data found in /Data/BC

This is the breast cancer data set
'''

import os
import glob

BC_DATA_DIR = os.getcwd() + "/Data/HybridSets/BC"

BC_EXPRESSION_DIR = BC_DATA_DIR + "/ExpressionProfiles"
BC_CLINICAL_DATA_FILE = BC_DATA_DIR + "/ClinicalData/ClinicalData.txt"

'''
Parameters:
file    =   which is incorrectly formatted
error   =   a custom error message
'''
class FileFormatError(Exception):
    def __init__(self, file, error):
        self.file = file
        self.error = error

    def __str__(self):
        return "File " + str(self.file) + "\n" + str(self.error)

'''
Wrapper for expression data, see BC_DATA_DIR/ExpressionProfiles for more details
There is a one-one correspondance with fields here and fields in the initial dataset
Don't really need to know what they are
'''
class expression_profile:
    def __init__(self, sample_num, substance, gene, log_ratio, log_error, p_value, intensity, flag):
        self.sample_num = sample_num
        self.substance = substance
        self.gene = gene
        self.log_ratio = log_ratio
        self.log_error = log_error
        self.p_value = p_value
        self.intensity = intensity
        self.flag = flag

    def __str__(self):
        st     =   "\tSAMPLE\t" + str(self.sample_num) + "\n"
        st     +=  "\tSUBSTANCE\t" + str(self.substance) + "\n"
        st     +=  "\tGENE\t" + str(self.gene) + "\n"
        st     +=  "\tLOG RATIO\t" + str(self.log_ratio) + "\n"
        st     +=  "\tLOG ERROR\t" + str(self.log_error) + "\n"
        st     +=  "\tp-Value\t" + str(self.p_value) + "\n"
        st     +=  "\tINTENSITY\t" + str(self.intensity) + "\n"
        st     +=  "\tFLAG\t" + str(self.flag) + "\n"
        return st

'''
As expression_profile, but for clinical_data
'''
class clinical_data:
    def __init__(self, sample_num, first_series, posnodes, event_meta,
                 event_death, time_survival, time_recur, time_meta, esr1, nih, st_gallen, conserv,
                 c1_from_data, c1_cross_valid, c1_used):
        self.sample_num = sample_num
        self.first_series = first_series
        self.posnodes = posnodes
        self.event_meta = event_meta
        self.event_death = event_death
        self.time_survival = time_survival
        self.time_recur = time_recur
        self.time_meta = time_meta
        self.esr1 = esr1
        self.nih = nih
        self.st_gallen = st_gallen
        self.conserv = conserv
        self.c1_from_data = c1_from_data
        self.c1_cross_valid = c1_cross_valid
        self.c1_used = c1_used


class sample:
    def __init__(self, profiles, sample_num):
        self.profiles = profiles
        self.sample_num = sample_num

'''
Given a file containing expression profiles, a list of expression_profiles based on the data
See /Data/HybridSets/BC/ExpressionProfiles for the format of data
'''
def readExpressionProfile(file_name):
    file = open(file_name, 'r')

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

        #some of the data has a random newline ¯\_(ツ)_/¯
        if str(row) == '':
            break

        #each row contains multiple samples for each gene45
        fields = row.split("\t")

        substance = fields[0]
        gene = fields[1]

        #ToDo: remove magic numbers
        #start from two, which is the index of the first column with actual information
        #5 fields, so offset by 5 with each id
        for f in range(0, len(IDs)):
            sample_num = IDs[f]
            log_ratio = fields[2 + 5 * f + 0]
            log_error = fields[2 + 5 * f + 1]
            p_value = fields[2 + 5 * f + 2]
            intensity = fields[2 + 5 * f + 3]
            flag = fields[2 + 5 * f + 4]

            profile = expression_profile(sample_num=sample_num,
                                         substance=substance,
                                         gene=gene,
                                         log_ratio=log_ratio,
                                         log_error=log_error,
                                         p_value=p_value,
                                         intensity=intensity,
                                         flag=flag)
            expression_profiles.append(profile)

    file.close()
    return expression_profiles

'''
Assumes all text files in dir are .txt files formatted as per what BC_DATA_DIR/ExpressionProfiles/README says

dir     =       the directory containing the .txt files containing expression profiles

returns a formatted list of expression profiles
'''
def getExpressionProfiles(dir):
    files = [filename for filename in glob.glob(os.path.join(dir, "*.txt"))]
    expression_profiles = []
    for file in files:
        expression_profiles += readExpressionProfile(file)

    return expression_profiles

'''
Given a file containing clinical data, return a list of expression_profiles based on the data
See /Data/HybridSets/BC/ClinicalData for the format of data
'''
def getClinicalData(file_name):
    file = open(file_name, 'r')

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

        data = clinical_data(sample_num, first_series, posnodes, event_meta,
                             event_death, time_survival, time_recur, time_meta,
                             esr1, nih, st_gallen, conserv,
                             c1_from_data, c1_cross_valid, c1_used)
        clinical_datas.append(data)

    return clinical_datas