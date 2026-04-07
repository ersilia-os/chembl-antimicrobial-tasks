import os
import sys
import zipfile
import pandas as pd
import numpy as np

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from model_utils import load_all_gz_csvs_from_zip


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
    datasets_path = os.path.join(base_path, "12_datasets")
    final_datasets_17 = os.path.join(base_path, "17_final_datasets.csv")
    output_zip = os.path.join(base_path, "19_final_datasets.zip")
    metadata_csv = os.path.join(base_path, "19_final_datasets_metadata.csv")

    if not os.path.exists(final_datasets_17):
        raise FileNotFoundError(f"Required input file not found: {final_datasets_17}")

    # Selected datasets are determined by step 17's 'selected' column
    print(f"Reading selected datasets from {final_datasets_17}")
    final_17 = pd.read_csv(final_datasets_17)
    selected_rows = final_17[final_17['selected'] == True].reset_index(drop=True)
    print(f"Found {len(selected_rows)} selected datasets")

    if len(selected_rows) == 0:
        print("No selected datasets found. Nothing to process.")
        return

    datasets_to_export = []
    for _, row in selected_rows.iterrows():
        is_merged = row['name'].startswith('M_')
        datasets_to_export.append({
            'name': row['name'],
            'label': row['label'],
            'activity_type': row.get('activity_type', ''),
            'unit': row.get('unit', ''),
            'target_type': row.get('target_type', ''),
            'cutoff': row.get('cutoff', ''),
            'auroc': row.get('auroc', np.nan),
            'cpds': row.get('cpds', 0),
            'positives': row.get('positives', 0),
            'source': 'merged' if is_merged else 'individual',
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
                simplified_df = (simplified_df
                    .sort_values("bin", ascending=False)
                    .drop_duplicates(subset="smiles", keep="first")
                    .reset_index(drop=True))
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