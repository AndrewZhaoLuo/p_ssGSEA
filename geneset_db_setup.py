'''
This file contains how to import gene_sets into database. Note, there is no support for seperate genesets for each
dataset. If this is needed in the future, send an email
'''

import os
import sqlite3

GENE_SET_DIR = "./Data/AppCache/"
GENE_SET_DB = GENE_SET_DIR + "GeneSets.db"

def create_GeneSet_schema_sets(cursor):
    cursor.execute('''CREATE TABLE GeneSet_URL
        (GeneSet text, URL text)''')

    cursor.execute('''CREATE TABLE GeneSet_Genes
        (GeneSet text, Gene text)''')

def load_gene_sets(cursor, gene_sets):
    if not os.path.exists(GENE_SET_DIR):
        os.makedirs(GENE_SET_DIR)

    for set in gene_sets:
        set_name = set.set_name
        set_url = set.url
        cursor.execute("INSERT INTO GeneSet_URL VALUES(?,?)",
                       (set_name, set_url))

        genes = set.genes
        for gene in genes:
            cursor.execute("INSERT INTO GeneSet_Genes VALUES(?,?)",
                           (set_name, gene))

def create_GeneSet_db(sets):
    '''
    Creates a database with gene set data with the given gene_sets

    :param sets: a list of gene_sets
    :type sets: list
    '''

    if not os.path.exists(GENE_SET_DIR):
        os.makedirs(GENE_SET_DIR)

    cursor = sqlite3.connect(GENE_SET_DB)

    create_GeneSet_schema_sets(cursor)
    load_gene_sets(cursor, sets)

    cursor.commit()
    cursor.close()


if __name__ == "__main__":
    import process_GeneSets

    #read some sets
    sets = process_GeneSets.readAllGeneSets(process_GeneSets.GENE_SETS_DIR)

    #load! ez pz!!!
    create_GeneSet_db(sets)