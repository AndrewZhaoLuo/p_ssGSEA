'''
This file contains general purpose for importing expression data onto a database
'''

import sqlite3
import os
from process_BC_data import readAllExpressionProfiles
from process_BC_data import BC_EXPRESSION_DIR

#general
DB_BASE_DIR = lambda data_set: os.getcwd() + "/Data/AppCache/" + data_set + "/"
DB_GENE_EXPRESSION = lambda data_set: data_set + "_GeneExpression.db"
TABLE_GENE_EXPRESSION = lambda  data_set: data_set + "_GeneExpression"

#Schema
def create_schema_expression(cursor, other_fields, data_set):

    query = "CREATE TABLE " + data_set + "_GeneExpression(Sample int, Gene text, Intensity real)"
    cursor.execute(query)

    query = lambda field: "ALTER TABLE " + TABLE_GENE_EXPRESSION(data_set) +" ADD COLUMN " + field + " text"
    for field in other_fields:
        cursor.execute(query(field))

def load_db_expression(cursor, profiles, data_set):

    for profile in profiles:
        other_fields = profile.other_fields

        id = profile.id
        gene = profile.gene
        intensity = profile.intensity

        command_tuple = (id, gene, intensity)
        for field in other_fields:
            command_tuple += (other_fields[field],)

        query = "INSERT INTO " + TABLE_GENE_EXPRESSION(data_set) + " VALUES(?,?,?,?,?,?,?,?)"
        cursor.execute(query, command_tuple)

def create_expression_db(data_set, profiles):
    '''
    Creates a database with expression data with the given data_set name from the list of given profiles.

    :param data_set: the name of the data_set, creates own folder to do all your cool analysis in!!
    :type data_set: str

    :param profiles: a list of expression_profiles
    :type profiles: list
    '''

    #Create directory
    if not os.path.exists(DB_BASE_DIR(data_set)):
        os.makedirs(DB_BASE_DIR(data_set))

    cursor = sqlite3.connect(DB_BASE_DIR(data_set) + DB_GENE_EXPRESSION(data_set))

    create_schema_expression(cursor, profiles[0].other_fields, data_set)
    load_db_expression(cursor, profiles, data_set)

    cursor.commit()
    cursor.close()


if __name__ == "__main__":
    #example use:
    #load up data formatted into expression_profiles
    profiles = readAllExpressionProfiles(BC_EXPRESSION_DIR)

    #create label and go!
    label = "BC2"
    create_expression_db(label, profiles)