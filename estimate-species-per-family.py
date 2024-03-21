import csv
import argparse
from collections import defaultdict, Counter
from sourmash.tax.tax_utils import BaseLineageInfo


ICTV_ranks = ('realm','subrealm','kingdom','subkingdom','phylum','subphylum',
              'class','subclass','order','suborder','family','subfamily',
              'genus','subgenus','species','name')

# Function to count unique names in each family, and unique species and genera
def parse_lineage_file(csv_file_path):
    # first, make family dictionary
    rankLineageD = defaultdict(lambda: Counter())
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            # full_lineage = row['lineages'] + ';' + row['name']
            species_lineage = row['lineage']
            lineage_info = BaseLineageInfo(ICTV_ranks, lineage_str=species_lineage)
            family_lineage = lineage_info.pop_to_rank('family')
            # species_name = lineage_info.name_at_rank('species')
            rankLineageD[family_lineage][species_lineage] += 1 # keep track of how many times we saw each species
    
    return rankLineageD


# Main function to process the file and output results
def main(args):
    rankLineageDict = parse_lineage_file(args.lineage_csv)
    
    with open(args.output_csv, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['family_lineage', 'n_unique_species'])
        # writer.writerow(['family_lineage', 'unique_species_counts', 'species_with_counts'])

        for family_lin, species_counter in rankLineageDict.items():
            n_unique_species = len(species_counter)
            # species_counts_str = ";".join([f"{species}: {count}" for species, count in species_counter.items()])
            writer.writerow([family_lin.display_lineage(), n_unique_species])
            # writer.writerow([family_lin.display_lineage(), unique_species_count, species_counts_str])

    with open(args.species_counts, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['species_lineage', 'count'])

        for family_lin, species_counter in rankLineageDict.items():
            for species_lin, count in species_counter.items():
                writer.writerow([species_lin, count])


# Command line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process lineage information from a CSV file.")
    parser.add_argument("lineage_csv", help="CSV file with lineage information")
    parser.add_argument("-o", "--output-csv", help="Output CSV file to save the results")
    parser.add_argument('-s', '--species-counts', help='output counts per species')
    args = parser.parse_args()
    main(args)