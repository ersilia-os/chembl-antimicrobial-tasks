import os
import sys
import zipfile
import pandas as pd
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from model_utils import load_all_gz_csvs_from_zip


def _norm_key(key):
    """Normalize NaN unit to None so the tuple is hashable and equality works."""
    a, at, u = key
    return (a, at, None if (not isinstance(u, str) and pd.isna(u)) else u)


def _parse_assay_key(s):
    """Parse assay key string into (assay_id, activity_type, unit) tuple."""
    assay_id, activity_type, unit = s.split("|")
    return (assay_id, activity_type, np.nan if unit == "" else unit)


def main(pathogen_code: str) -> None:
    """
    Create final datasets using selected assays from table 18 with original names from tables 16 and 17.

    Parameters
    ----------
    pathogen_code : str
        Pathogen code (e.g., 'mtuberculosis')
    """
    print(f"Preparing final datasets for {pathogen_code}...")

    # Define paths
    base_path = os.path.join(os.path.dirname(__file__), "..", "output", pathogen_code)
    datasets_path = os.path.join(base_path, "datasets")
    assay_master_csv = os.path.join(base_path, "18_assays_master.csv")
    final_datasets_17 = os.path.join(base_path, "17_final_datasets.csv")
    output_zip = os.path.join(base_path, "19_final_datasets.zip")
    metadata_csv = os.path.join(base_path, "19_final_datasets_metadata.csv")

    # Check required files
    for required_file in [assay_master_csv, final_datasets_17]:
        if not os.path.exists(required_file):
            raise FileNotFoundError(f"Required input file not found: {required_file}")

    # Read table 18 to get all selected assays
    print(f"Reading selected datasets from {assay_master_csv}")
    assay_master = pd.read_csv(assay_master_csv)
    selected_assays = assay_master[assay_master['selected'] == True].copy()
    print(f"Found {len(selected_assays)} selected assays")

    if len(selected_assays) == 0:
        print("No selected datasets found. Nothing to process.")
        return

    # Create a set of selected assay keys for fast lookup
    selected_keys = set()
    for _, row in selected_assays.iterrows():
        key = _norm_key((row['assay_id'], row['activity_type'], row['unit']))
        selected_keys.add(key)

    print(f"Created lookup set with {len(selected_keys)} selected assay keys")

    # Get datasets to export from table 17 (contains all final selected datasets)
    datasets_to_export = []

    # Get all datasets from table 17 (both individual and merged)
    print(f"Reading all final dataset names from {final_datasets_17}")
    final_17 = pd.read_csv(final_datasets_17)

    for _, row in final_17.iterrows():
        # Check if any assays in this dataset are selected
        dataset_selected = False
        for assay_key_str in row['assay_keys'].split(';'):
            assay_key = _norm_key(_parse_assay_key(assay_key_str))
            if assay_key in selected_keys:
                dataset_selected = True
                break

        if dataset_selected:
            # Determine if this is a merged dataset based on name
            is_merged = row['name'].startswith('M_')
            source = 'merged' if is_merged else 'individual'

            datasets_to_export.append({
                'name': row['name'],
                'label': row['label'],
                'activity_type': row.get('activity_type', 'MULTIPLE' if is_merged else ''),
                'unit': row.get('unit', 'MULTIPLE' if is_merged else ''),
                'target_type': row.get('target_type', ''),
                'cutoff': row.get('cutoff', ''),
                'auroc': row.get('auroc', np.nan),
                'cpds': row.get('cpds', 0),
                'positives': row.get('positives', 0),
                'source': source
            })

    print(f"Found {len(datasets_to_export)} datasets to export")

    if len(datasets_to_export) == 0:
        print("No datasets found to export.")
        return

    # Load all available datasets
    all_datasets = {}

    # Load individual datasets from ZIP files
    zip_files = {
        'quantitative': os.path.join(datasets_path, "datasets_qt.zip"),
        'qualitative': os.path.join(datasets_path, "datasets_ql.zip"),
        'mixed': os.path.join(datasets_path, "datasets_mx.zip")
    }

    for dataset_type, zip_path in zip_files.items():
        if os.path.exists(zip_path):
            print(f"Loading {dataset_type} datasets from {zip_path}")
            datasets_dict = load_all_gz_csvs_from_zip(zip_path)
            all_datasets.update(datasets_dict)

    # Load merged datasets from M/ directory
    merged_dir = os.path.join(datasets_path, "M")
    if os.path.exists(merged_dir):
        print(f"Loading merged datasets from {merged_dir}")
        for filename in os.listdir(merged_dir):
            if filename.endswith('.csv.gz'):
                filepath = os.path.join(merged_dir, filename)
                try:
                    dataset = pd.read_csv(filepath)
                    all_datasets[filename] = dataset
                except Exception as e:
                    print(f"Warning: Could not load {filename}: {e}")

    print(f"Loaded {len(all_datasets)} total datasets")

    # Export selected datasets using original names
    processed_count = 0
    metadata_records = []

    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as out_zip:
        for dataset_info in datasets_to_export:
            original_name = dataset_info['name']
            dataset_key = f"{original_name}.csv.gz"

            if dataset_key in all_datasets:
                dataset = all_datasets[dataset_key]

                # Verify required columns
                if 'smiles' not in dataset.columns or 'bin' not in dataset.columns:
                    print(f"Warning: Dataset {original_name} missing required columns (smiles, bin)")
                    continue

                # Export with original name (just remove .gz and keep .csv)
                simplified_df = dataset[['smiles', 'bin']].copy()
                csv_content = simplified_df.to_csv(index=False)
                out_zip.writestr(f"{original_name}.csv", csv_content)

                # Build metadata record
                metadata_records.append({
                    'original_name': original_name,
                    'activity_type': dataset_info['activity_type'],
                    'unit': dataset_info['unit'],
                    'target_type': dataset_info['target_type'],
                    'cutoff': dataset_info['cutoff'],
                    'auroc': dataset_info['auroc'],
                    'cpds': len(dataset),
                    'positives': dataset['bin'].sum(),
                    'label': dataset_info['label'],
                    'source': dataset_info['source']
                })

                processed_count += 1
                print(f"Exported: {original_name}")

            else:
                print(f"Warning: Dataset file not found: {dataset_key}")

    # Save metadata CSV
    metadata_df = pd.DataFrame(metadata_records)
    metadata_df.to_csv(metadata_csv, index=False)

    print(f"Successfully exported {processed_count} datasets")
    print(f"Created final datasets ZIP: {output_zip}")
    print(f"Saved metadata to: {metadata_csv}")
    print("All datasets exported with original names (no renaming)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/19_prepare_final_datasets.py <pathogen_code>")
        print("Example: python scripts/19_prepare_final_datasets.py mtuberculosis")
        sys.exit(1)

    pathogen_code = sys.argv[1]
    try:
        main(pathogen_code)
        print(f"✓ Successfully prepared final datasets for {pathogen_code}")
    except Exception as e:
        print(f"Error processing {pathogen_code}: {e}")
        raise