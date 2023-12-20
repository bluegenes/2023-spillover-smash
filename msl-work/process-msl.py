import argparse
import os
import csv
import yaml
import pandas as pd
import requests
from io import BytesIO
import sourmash
#from sourmash.tax_utils import

def download_file(url, name, year, out_dir, force=False):
    # Create a filename using the msl name and year
    filename = f"{name}_{year}.xlsx"
    filepath = os.path.join(out_dir, filename)

    # Check if the file already exists
    if not os.path.exists(filepath):
        print(f"Downloading {filename}...")
        response = requests.get(url)
        response.raise_for_status()

        # Write the content to the file
        with open(filepath, 'wb') as f:
            f.write(response.content)
    else:
        print(f"File {filename} already exists. Skipping download.")

    return filepath


def process_msl_tab(df):
    # Current ICTV ranks (consider past within the 15-rank system)
    columns = ["Realm", "Subrealm", "Kingdom", "Subkingdom", "Phylum", "Subphylum",
               "Class", "Subclass", "Order", "Suborder", "Family", "Subfamily",
               "Genus", "Subgenus", "Species"]

    # Process each row to create a ';'-separated list
    processed_data = []
    for _, row in df.iterrows():
        # Create a list of values for each column, using an empty string for missing or null values
        values = [str(row[col]) if col in df and pd.notnull(row[col]) else '' for col in columns]
        concatenated_values = ';'.join(values)
        processed_data.append(concatenated_values)
    return processed_data


def find_matching_columns(df, keyword):
    # Find columns matching the keyword
    return [col for col in df.columns if keyword in col]

def process_changelog_tab(df):
    # Find columns that contain 'Old Name' and 'New Name'
    old_name_cols = find_matching_columns(df, 'Old Name')
    new_name_cols = find_matching_columns(df, 'New Name')

    # Old Name:: New Name
    name_mapping = {}
    for new_col, old_col in zip(new_name_cols, old_name_cols):
        # Skip if either column is missing
        if new_col not in df or old_col not in df:
            continue

        # Update the dictionary with new names as keys and old names as values
        name_mapping.update(df.set_index(old_col)[new_col].dropna().to_dict())

    return name_mapping

#def export_name_mapping_to_csv(name_mapping, output_path):
#    # Convert the name mapping to a DataFrame and write to a CSV file
#    name_mapping_df = pd.DataFrame(list(name_mapping.items()), columns=['old-lineage', 'new-lineage'])
#    name_mapping_df.to_csv(output_path, index=False)

def export_name_mapping_to_csv(name_mapping, output_path):
    with open(output_path, mode='w', newline='', encoding='utf-8') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write the header
        csvwriter.writerow(['old-lineage', 'new-lineage'])

        # Write the data
        for old_lineage, new_lineage in name_mapping.items():
            csvwriter.writerow([old_lineage, new_lineage])

def write_lineages_to_csv(lineages_dict, output_csv_path):
    with open(output_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write the header
        csvwriter.writerow(['msl_name', 'lineage'])

        # Write the data
        for msl_name, lineages in lineages_dict.items():
            for lineage in lineages:
                csvwriter.writerow([msl_name, lineage])

def main(args):
    with open(args.config, 'r') as file:
        config = yaml.safe_load(file)

    mslD = {}
    change_mapping = {}
    for file_config in config.get('msl_files', []):
        msl_name = file_config['name']
        print(f"Processing file: {file_config['name']} ({file_config['year']}) from {file_config['url']}")
        filepath = download_file(file_config['url'], msl_name, file_config['year'], args.download_directory, args.force)
        excel_data = pd.ExcelFile(filepath)

        # Process MSL tab
        msl_df = pd.read_excel(excel_data, sheet_name=file_config['msl'])
        msl_lineages = process_msl_tab(msl_df)
        mslD[msl_name] = msl_lineages

        # Process ChangeLog tab, if exists
        ch_tab = file_config.get('changelog')
        if ch_tab and ch_tab in excel_data.sheet_names:
            print(f"Read changeLog data from MSL '{msl_name}'")
            changelog_df = pd.read_excel(excel_data, sheet_name=file_config['changelog'])
            # likely only need to look at species level -- all other changes should be incorporated in species lineage
            changelog_spp_only = changelog_df[changelog_df['Rank'] == 'species']
            these_changes = process_changelog_tab(changelog_spp_only)
            change_mapping.update(these_changes)
            #print(f"ChangeLog data (dict): {changelog_dict}")

    # write csvs
    export_name_mapping_to_csv(change_mapping, args.output_msl_changes)
    write_lineages_to_csv(mslD, args.output_msl)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Excel files based on a YAML configuration.")
    parser.add_argument('config', help="Path to the YAML configuration file")
    parser.add_argument('-d', '--download-directory', help='path to output all downloaded MSL files')
    parser.add_argument('-o', '--output-msl', help='path to output aggregated msl csv')
    parser.add_argument('-c', '--output-msl-changes', help='path to output aggregated lineage changes csv')
    parser.add_argument('-f', '--force', help='download files even if they already exist (force overwrite)')
    args = parser.parse_args()
    main(args)

