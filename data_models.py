'''
This file contains the models used for structuring data in gene expression data
'''

'''
Contains data representing a single output of expression of a single gene
Parameters:
id          =       an identifier unique to the person this sample was taken from
gene        =       the name of the gene as a string, suggested to use name provided by the
                        HUGO Gene Nomenclature Committee
intensity   =       intensity of the gene expressed as a real number, suggested units be FPKM
other_fields=       map of additional labelled information mapped ot their values. key should be the name of the field
'''
class expression_profile:
    def __init__(self, id, gene, intensity, other_fields):
        self.id = id
        self.gene = gene
        self.intensity = intensity
        self.other_fields = other_fields
'''
Contains clinical data (ie phenotypes) for patients who contributed genetic data for expression study

Parameters:
id          =       an identifier unique to the person this sample was taken from
other_fields=       map of additional labelled information mapped ot their values. key should be the name of the fields
'''
class clinical_data:
    def __init__(self, id, other_fields):
        self.id = id
        self.other_fields = other_fields

'''
A collection of genes pre-grouped due to being on similiar pathway, etc.

set_name    =       name of the pathway assoc. with this gene set
url         =       broad institute link to pathway info
genes       =       a list-like collection of strings representing genes in this set
'''
class gene_set:
    def __init__(self, set_name, url, genes):
        self.set_name = set_name
        self.url = url
        self.genes = genes

'''
A patient with a collection of all his/her gene expression profiles

profiles    =       a map of expression_profiles for this sample subject where key = gene name, value = expression_profile
sample_num  =       the sample's unique number identifier
'''
class sample:
    def __init__(self, profiles, sample_num):
        self.profiles = profiles
        self.sample_num = sample_num