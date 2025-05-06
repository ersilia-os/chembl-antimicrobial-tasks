#!/bin/bash
set -e  # Exit immediately if any command fails

# Define the list of pathogens
PATHOGENS=("abaumannii" "calbicans" "campylobacter" "ecoli" "efaecium" "enterobacter" "hpylori" "kpneumoniae" "mtuberculosis" "ngonorrhoeae" "paeruginosa" "pfalciparum" "saureus" "smansoni" "spneumoniae")
PATHOGENS=("abaumannii" "kpneumoniae")


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

    # Conditionally run 013a or 013b based on task type
    if [ "$TASK_TYPE" == "organism" ]; then
        python $SCRIPT_DIR/../src/013a_binarize_fetched_pathogen_data_ORG.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
    elif [ "$TASK_TYPE" == "protein" ]; then
        python $SCRIPT_DIR/../src/013b_binarize_fetched_pathogen_data_SP.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
    fi

    python $SCRIPT_DIR/../src/014_datasets_modelability.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR --$TASK_TYPE  # Task type is organism or protein
    python $SCRIPT_DIR/../src/015_datasets_distinguishability.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR --$TASK_TYPE  # Task type is organism or protein
    python $SCRIPT_DIR/../src/016_select_tasks_MOD_DIS.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR
    python $SCRIPT_DIR/../src/017_select_tasks_RED.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR --$TASK_TYPE  # Task type is organism or protein
    python $SCRIPT_DIR/../src/018_calculate_correlations.py --pathogen_code $PATHOGEN_CODE --output_dir $OUTPUT_DIR --$TASK_TYPE  # Task type is organism or protein


done
