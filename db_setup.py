'''
This is a simple series of script which sets up and retrieves information from the database

IN GENERAL THIS MESSES WITH THE DB SCHEMA, DONT MESS WITH UNLESS YOU KNOW WHAT YOU ARE DOING!!!
'''

import process_BC_data as BC

import timeit

'''
DATABASE CODE
'''
#Schema
def create_BC_schema_expression(cursor):
    cursor.execute('''CREATE TABLE BC_GeneExpression
        (Sample int, Substance text, Gene text, Log_Ratio real,
         Log_Ratio_Error real, p_Value real, Intensity real, Flag int)''')

def create_BC_schema_clinical(cursor):
    cursor.execute('''CREATE TABLE BC_ClinicalData
        (Sample int, FirstSeriesID int, Posnodes text, EventMeta int,
        EventDeath int, TimeSurvival real, TimeRecur real, TimeMeta real, ESR1 int,
        NIH int, StGallen int, Conserv int, C1FromData real, C1CrossValid real,
        C1Used real)''')

def create_BC_schema_sets(cursor):
    cursor.execute('''CREATE TABLE BC_GeneSet_URL
        (GeneSet text, URL text)''')

    cursor.execute('''CREATE TABLE BC_GeneSet_Genes
        (GeneSet text, Gene text)''')

def create_BC_schema_all(cursor):
    create_BC_schema_expression(cursor)
    create_BC_schema_clinical(cursor)
    create_BC_schema_sets(cursor)

#data loading
def load_BC_data_expression(cursor):
    gene_profiles = BC.getExpressionProfiles(BC.BC_EXPRESSION_DIR)
    for profile in gene_profiles:
        other_fields = profile.other_fields
        cursor.execute("INSERT INTO BC_GeneExpression VALUES(?,?,?,?,?,?,?,?)",
                       (profile.id, other_fields["substance"], profile.gene,
                        other_fields["log_ratio"], other_fields["log_error"], other_fields["p_value"],
                        profile.intensity, other_fields["flag"]))

def load_BC_data_clinical(cursor):
    clinical_profiles = BC.getClinicalData(BC.BC_CLINICAL_DATA_FILE)
    for profile in clinical_profiles:
        other_fields = profile.other_fields
        cursor.execute("INSERT INTO BC_ClinicalData VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                       (profile.id, other_fields["first_series"],other_fields["posnodes"], other_fields["event_meta"],
                        other_fields["event_death"], other_fields["time_survival"], other_fields["time_recur"], other_fields["time_meta"],
                        other_fields["esr1"], other_fields["nih"], other_fields["st_gallen"], other_fields["conserv"],
                        other_fields["c1_from_data"], other_fields["c1_cross_valid"], other_fields["c1_used"]))

def load_BC_data_sets(cursor):
    gene_sets = BC.getGeneSetData(BC.BC_GENE_SETS_FILE)
    for set in gene_sets:
        set_name = set.set_name
        set_url = set.url
        cursor.execute("INSERT INTO BC_GeneSet_URL VALUES(?,?)",
                       (set_name, set_url))

        genes = set.genes
        for gene in genes:
            cursor.execute("INSERT INTO BC_GeneSet_Genes VALUES(?,?)",
                           (set_name, gene))

def load_BC_data_all(cursor):
    load_BC_data_expression(cursor)
    load_BC_data_clinical(cursor)
    load_BC_data_sets(cursor)

'''
Rebuild BC database
'''
def rebuild_BC_db(cursor):
    print("Rebuilding BC database...")
    start = timeit.default_timer()

    #cursor.execute("DROP TABLE IF EXISTS BC_GeneExpression")
    #cursor.execute("DROP TABLE IF EXISTS BC_ClinicalData")
    #cursor.execute("DROP TABLE IF EXISTS BC_GeneSet_URL")
    #cursor.execute("DROP TABLE IF EXISTS BC_GeneSet_Genes")

    create_BC_schema_all(cursor)
    load_BC_data_all(cursor)

    end = timeit.default_timer()
    print("Finished building DB! Took " + str(end - start) + "s")
