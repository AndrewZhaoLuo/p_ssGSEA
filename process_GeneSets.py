'''
This file contains methods for importing and pre-processing the gene sets from MSigDB
'''

import os
import glob
from data_models import gene_set

GENE_SETS_DIR = os.getcwd() + "/Data/GeneSets/"

def readtGeneSetData(file_name):
    """
    Given a file and directory, reads a file containing gene set information assuming it is formatted as a .gmt

    :param file_name: the file to read gene information from
    :type file_name: str

    :returns: a list of gene_sets
    """
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
            #sanitize by removing / which can cause directory errors
            genes.append(set_info[i].replace("/", "-"))

        gene_sets.append(gene_set(set_name, set_url, genes))

    return gene_sets

def readAllGeneSets(dir):
    """
    Reads all files containing expression data in the given directory

    :param dir: the directory from where to read all gmt files containing expression data.
    :type dir: str

    :return: a list of gene_sets formed from all files read in the directory
    """
    files = [filename for filename in glob.glob(os.path.join(dir, "*.gmt"))]
    gene_sets = []
    for file in files:
        gene_sets += readtGeneSetData(file)

    return gene_sets