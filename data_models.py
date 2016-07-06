'''
This file contains the models used for structuring data in gene expression data
'''
class expression_profile:
    """
    Represents a single expression value of a single gene from a single sample.

    :param id: represents the unique identifier given to the person this sample was taken from
    :type id: int

    :param gene: the name of the gene whose expression is recorded in this profile
    :type gene: str

    :param intensity: the expression levels of the gene
    :type intensity: float

    :param other_fields: other fields associated with this profile. This takes the form of a map where the name
                         of each additional field is a key for the value of the field. The map should have strings
                         as keys going to objects as values
    :type other_fields: dict
    """
    def __init__(self, id, gene, intensity, other_fields):
        self.id = id
        self.gene = gene
        self.intensity = intensity
        self.other_fields = other_fields

class clinical_data:
    """
    Represents clinical data about each person where a gene sample was collected from

    :param id: represents the unique identifier given to the person this sample was taken from
    :type id: int

    :param other_fields: the fields associated with this person. This takes the form of a map where the name
                         of each additional field is a key for the value of the field. The map should have strings
                         as keys going to objects as values
    :type other_fields: dict
    """
    def __init__(self, id, other_fields):
        self.id = id
        self.other_fields = other_fields

class gene_set:
    """
    Represents an a priori gene set taken from the Broad institute's molecular signatures database

    :param set_name: the name of this set
    :type set_name: str

    :param url: the MSigDB webpage associated with this gene set
    :type url: str

    :param genes: a set of genes which are located in this set
    :type genes: set
    """

    def __init__(self, set_name, url, genes):
        self.set_name = set_name
        self.url = url
        self.genes = genes

class sample:
    """
    A grouping of gene expression profiles under the person who gave the gene sample

    :param profiles: A mapping of gene_names to expression_profiles
    :type profiles: dict

    :param id: the unique identifier for the person who gave samples which were used to derive the expression_profiles
    :type id: int
    """
    def __init__(self, profiles, id):
        self.profiles = profiles
        self.id = id