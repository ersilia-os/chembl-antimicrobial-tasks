import pandas as pd
import zipfile
import json
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))

# List of pathogens to process
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]
pathogens = ["Acinetobacter baumannii", "Mycobacterium tuberculosis", "Klebsiella pneumoniae"]

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

# For each pathogen
for pathogen in pathogens[1:2]:

    # Get pathogen code
    pathogen_code = get_pathogen_code(pathogen)

    # Define path to parameters
    PATH_TO_PARAMETERS = os.path.join(OUTPUT, pathogen_code, "assay_parameters.zip")

    # Load assays info
    ASSAYS_CLEANED = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_cleaned.csv"))
    ASSAYS_CLUSTERS = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_clusters.csv"))
    ASSAYS_DATA_RANGES = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "assays_data.csv"))

    # Shared columns
    KEYS = ["assay_id", "activity_type", "unit"]

    # Columns to take from each table
    COLUMNS_CLEANED = ["assay_id", "assay_type", "assay_organism", "doc_chembl_id", "target_type", "target_chembl_id", "target_organism", "activity_type", 
                    "unit", "canonical_unit", "activities", "nan_values", "cpds", "direction", "activity_comment_counts", "standard_text_count"]
    COLUMNS_CLUSTERS = ['clusters_0.3', 'clusters_0.6', 'clusters_0.85']
    COLUMNS_DATA_RANGES = ["equal", 'higher', 'lower', "dataset_type", "expert_cutoff", "pos_qt", "ratio_qt", "cpds_qt", "min_", "p1", "p25", "p50", "p75", 
                        "p99", "max_", "pos_ql", "ratio_ql", "cpds_ql"]
    
    def assert_unique_on_keys(df, keys, name):
        dups = df.duplicated(keys, keep=False)
        if dups.any():
            examples = df.loc[dups, keys].drop_duplicates().head(10)
            raise ValueError(
                f"{name} is NOT unique on keys {keys} (found {dups.sum()} duplicate rows).\n"
                f"Example duplicate keys (first 10):\n{examples.to_string(index=False)}")
        
    # Check that [assay, activity_type, unit] items are unique
    assert_unique_on_keys(ASSAYS_CLEANED, KEYS, "ASSAYS_CLEANED")
    assert_unique_on_keys(ASSAYS_CLUSTERS, KEYS, "ASSAYS_CLUSTERS")
    assert_unique_on_keys(ASSAYS_DATA_RANGES, KEYS, "ASSAYS_DATA_RANGES")

    ORGNISM_CURATED, TARGET_TYPE_CURATED, STRAIN, ATCC_ID, MUTATIONS, KDR, MEDIA = [], [], [], [], [], [], []

     # Inside zip file
    with zipfile.ZipFile(PATH_TO_PARAMETERS) as z:

        # Iterating over assays
        for assay_id, activity_type, unit in ASSAYS_CLEANED[['assay_id', 'activity_type', 'unit']].values:

            # Prepare filename
            filename = "_".join([str(assay_id), str(activity_type), str(unit), 'parameters']) + ".json"
            
            # Read JSON file inside zip
            with z.open(filename) as file:
                par = json.load(file)

            # Store results
            ORGNISM_CURATED.append(par['organism'])
            TARGET_TYPE_CURATED.append(par['target_type_curated'])
            STRAIN.append(par['strain'])
            ATCC_ID.append(par['atcc_id'])
            MUTATIONS.append(";".join(par['mutations']))
            KDR.append(";".join(par['known_drug_resistances']))
            MEDIA.append(par['media'])    

    # Complete table
    ASSAYS_CLEANED['organism_curated'] = ORGNISM_CURATED
    ASSAYS_CLEANED['target_type_curated'] = TARGET_TYPE_CURATED
    ASSAYS_CLEANED['strain'] = STRAIN
    ASSAYS_CLEANED['atcc_id'] = ATCC_ID
    ASSAYS_CLEANED['mutations'] = MUTATIONS
    ASSAYS_CLEANED['known_drug_resistances'] = KDR
    ASSAYS_CLEANED['media'] = MEDIA

    # Merge tables
    ASSAYS_CLEANED = ASSAYS_CLEANED.merge(ASSAYS_CLUSTERS[KEYS + COLUMNS_CLUSTERS],on=KEYS, how="left", validate="1:1")
    ASSAYS_CLEANED = ASSAYS_CLEANED.merge(ASSAYS_DATA_RANGES[KEYS + COLUMNS_DATA_RANGES],on=KEYS, how="left", validate="1:1")

    COLUMNS = ["assay_id", 
           "assay_type", 
           "assay_organism",
           "organism_curated",
           "doc_chembl_id",
           "target_type", 
           "target_type_curated",
           "target_chembl_id",
           "target_organism",
           "strain",
           "atcc_id",
           "mutations",
           "known_drug_resistances",
           "media", 
           "activity_type", 
           "unit",
           "canonical_unit", 
           "activities", 
           "nan_values", 
           "cpds",	
           "direction", 
           "activity_comment_counts", 
           "standard_text_count",
           "equal",
           'higher',
           'lower',
           "dataset_type",
           "expert_cutoff",
           "pos_qt",
           "ratio_qt",
           "cpds_qt",
           "min_",
           "p1",
           "p25",
           "p50",
           "p75", 
           "p99", 
           "max_", 
           "pos_ql", 
           "ratio_ql", 
           "cpds_ql",
           "clusters_0.3",
           "clusters_0.6",
           "clusters_0.85"]
    
    # Sort columns
    ASSAYS_CLEANED = ASSAYS_CLEANED[COLUMNS]

    # Save as CSV
    ASSAYS_CLEANED.to_csv(os.path.join(OUTPUT, pathogen_code, "assays_master.csv"), index=False)