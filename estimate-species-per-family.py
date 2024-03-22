import csv
import argparse
from collections import defaultdict, Counter
from sourmash.tax.tax_utils import BaseLineageInfo, RankLineageInfo, ICTV_RANKS, ICTVRankLineageInfo, get_ident


# Function to count unique names in each family, and unique species and genera
def parse_fastmultigather_lineage_file(csv_file_path):
    annotationD = {} 
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ident = get_ident(row['name'])
            # TODO: check threshold is at least 3 x scaled, only keep annotation if passes threshold
            full_lineage = row['lineage'] + ';' + row['name']
            lineage_info = ICTVRankLineageInfo(lineage_str=full_lineage)
            annotationD[ident] = lineage_info
            # species_lineage = row['lineage'] # check -- is this species or virus level??
            # lineage_info = BaseLineageInfo(ICTV_RANKS, lineage_str=species_lineage)
            # lineage_info = ICTVRankLineageInfo(lineage_str=species_lineage)
            # species_name = lineage_info.name_at_rank('species')
    return annotationD

def parse_annotated_clusters(csv_file_path):
    cluster_annotD = {} 
    with open(csv_file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            lineage = row['lineage']
            idents = row['nodes'].split(';')
            for ident in idents:
                cluster_annotD[ident] = lineage
    return cluster_annotD


def count_per_family(annotationD, count_rank = 'species'):
    '''
    count the number of viruses (name) and count_rank (default species)
    per viral family
    '''
    virusLineageD = defaultdict(lambda: Counter())
    rankLineageD = defaultdict(lambda: Counter())

    for _, lineage_info in annotationD.items():
        family_lineage = lineage_info.pop_to_rank('family')
        rank_lineage = lineage_info.pop_to_rank(count_rank)
        virusLineageD[family_lineage][lineage_info] += 1 # keep track of how many times we saw each name
        rankLineageD[family_lineage][rank_lineage] += 1 # keep track of how many times we saw each RANK-level annot

    return rankLineageD


# Main function to process the file and output results
def main(args):
    # load all annotations
    dna_fmg = parse_fastmultigather_lineage_file(args.dna_lineages)
    prot_fmg = parse_fastmultigather_lineage_file(args.protein_lineages)
    dna_cluster = parse_annotated_clusters(args.dna_cluster_lineages)
    prot_cluster = parse_annotated_clusters(args.protein_cluster_lineages)

    # for each ident, check annotations from each source, combine as needed (LCA, pull back to genus/family?? for protein)
    annotD = {} # combined annotation dict
    with open(args.all_idents_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ident = get_ident(row['name'])

            dna_annot = dna_fmg.get(ident, '')
            dna_cluster_annot = dna_cluster.get(ident, '')
            prot_annot = prot_fmg.get(ident, '')
            prot_cluster_annot = prot_cluster(ident, '')

            # to do: combine annotations properly
            annotD[ident] = dna_annot
    
    virus_counts, rank_counts = count_per_family(annotD, args.rank)


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
    parser = argparse.ArgumentParser(description="Aggregate lineage information to final annotations.")
    parser.add_argument("all_idents_csv", help="CSV file with all identifiers for annotation")
    parser.add_argument("--dna-lineages", help="CSV file with dna-based lineage information")
    parser.add_argument("--protein-lineages", help="CSV file with protein-based lineage information")
    parser.add_argument("--dna-cluster-lineages", help="CSV file with dna-based annotated clusters")
    parser.add_argument("--protein-cluster-lineages", help="CSV file with protein-based annotated clusters")
    parser.add_argument("-o", "--output-csv", help="Output CSV file to save the results")
    parser.add_argument('-s', '--species-counts', help='output counts per species')
    parser.add_argument('-r', '--rank', default='species', help= 'rank to report counts')
    args = parser.parse_args()
    main(args)