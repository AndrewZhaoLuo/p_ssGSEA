'''
This is a simple series of script which

IN GENERAL THIS MESSES WITH THE DB SCHEMA, DONT MESS WITH UNLESS YOU KNOW WHAT YOU ARE DOING!!!
'''

import sqlite3
import process_BC_data as BC

import timeit

def create_BC_schema(cursor):
    cursor.execute('''CREATE TABLE BC_GeneExpression
        (Sample int, Substance text, Gene text, Log_Ratio real,
         Log_Ratio_Error real, p_Value real, Intensity real, Flag int)''')

    cursor.execute('''CREATE TABLE BC_ClinicalData
        (Sample int, FirstSeriesID int, Posnodes text, EventMeta int,
        EventDeath int, TimeSurvival real, TimeRecur real, TimeMeta real, ESR1 int,
        NIH int, StGallen int, Conserv int, C1FromData real, C1CrossValid real,
        C1Used real)''')

def load_BC_data(cursor):
    gene_profiles = BC.getExpressionProfiles(BC.BC_EXPRESSION_DIR)
    for profile in gene_profiles:
        cursor.execute("INSERT INTO BC_GeneExpression VALUES(?,?,?,?,?,?,?,?)",
                       (profile.sample_num, profile.substance, profile.gene,
                        profile.log_ratio, profile.log_error, profile.p_value,
                        profile.intensity, profile.flag))

    clinical_profiles = BC.getClinicalData(BC.BC_CLINICAL_DATA_FILE)
    for profile in clinical_profiles:
        cursor.execute("INSERT INTO BC_ClinicalData VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                       (profile.sample_num, profile.first_series, profile.posnodes, profile.event_meta,
                        profile.event_death, profile.time_survival, profile.time_recur, profile.time_meta,
                        profile.esr1, profile.nih, profile.st_gallen, profile.conserv,
                        profile.c1_from_data, profile.c1_cross_valid, profile.c1_used))

'''
Tears down and rebuild BC database (PLEASE DON'T CALL UNLESS YOU KNOW WHAT YOU ARE DOING!!!!!)
'''
def rebuild_BC_db(cursor):
    print("Rebuilding BC database...")
    start = timeit.default_timer()

    cursor.execute("DROP TABLE IF EXISTS BC_GeneExpression")
    cursor.execute("DROP TABLE IF EXISTS BC_ClinicalData")

    create_BC_schema(cursor)
    load_BC_data(cursor)

    end = timeit.default_timer()
    print("Finished building DB! Took " + str(end - start) + "s")

if __name__ == "__main__":
    connection = sqlite3.connect("GeneExpression.db")

    cursor = connection.cursor()
    #rebuild_BC_db(cursor)

    connection.commit()

