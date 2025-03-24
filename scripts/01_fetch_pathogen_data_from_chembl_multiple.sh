#!/bin/bash

# Define the list of pathogens
PATHOGENS=("abaumannii" "calbicans" "campylobacter" "ecoli" "efaecium" "enterobacter" "hpylori" "kpneumoniae" "mtuberculosis" "ngonorrhoeae" "paeruginosa" "pfalciparum" "saureus" "smansoni" "spneumoniae")

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# Function to display usage
usage() {
    echo "Usage: $0 --output_dir <directory>"
    echo "  --output_dir      Specify the output directory (absolute path)"
    exit 1
}

# Check for arguments
if [ $# -ne 2 ]; then
    echo "Error: Incorrect number of arguments."
    usage
fi

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --output_dir)
            OUTPUT_DIR=$2
            shift 2
            ;;
        *)
            echo "Error: Unknown argument $1"
            usage
            ;;
    esac
done

# Validate output directory
if [ -z "$OUTPUT_DIR" ]; then
    echo "Error: --output_dir is required."
    usage
fi
# Loop through each pathogen and run the scripts
for PATHOGEN_CODE in "${PATHOGENS[@]}"; do
    echo "Fetching data for pathogen $PATHOGEN_CODE and storing in $OUTPUT_DIR"
    
    python $SCRIPT_DIR/../src/011_pathogen_getter.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR 
    python $SCRIPT_DIR/../src/012_clean_fetched_pathogen_data.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
    python $SCRIPT_DIR/../src/013_binarize_fetched_pathogen_data.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
    python $SCRIPT_DIR/../src/014_datasets_modelability.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
    python $SCRIPT_DIR/../src/015_select_tasks.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
    python $SCRIPT_DIR/../src/016_wrapup_tasks_and_clean_output_folder.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR

done
