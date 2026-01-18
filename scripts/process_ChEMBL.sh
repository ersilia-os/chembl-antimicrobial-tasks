SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# python $SCRIPT_DIR/../src/00_export_chembl_activities.py
python $SCRIPT_DIR/../src/01_get_compound_info.py
# python $SCRIPT_DIR/../src/02_standardize_compounds.py
