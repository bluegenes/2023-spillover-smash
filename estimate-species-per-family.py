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
    lineageCountsD = defaultdict(lambda: {'species': Counter(), 'name': Counter()})

    for _, lineage_info in annotationD.items():
        family_lineage = lineage_info.pop_to_rank('family')
        rank_lineage = lineage_info.pop_to_rank(count_rank)
        lineageCountsD[family_lineage][count_rank][rank_lineage] += 1
        lineageCountsD[family_lineage]['virus'][lineage_info] += 1

    return lineageCountsD


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

            # to do: combine annotations properly (LCA as needed)
            annotD[ident] = dna_annot
    
    # count n unique virus, spp per family
    countsD = count_per_family(annotD, args.rank)
    with open(args.output_csv, 'w', newline='') as outfile:
        w = csv.writer(outfile)
        w.writerow(['family_lineage', 'n_species_counts', 'n_virus_counts'])
        for family_lineage, counts_dict in countsD.items():
            n_at_rank = len(counts_dict[args.rank])
            n_virus = len(counts_dict['name'])
            w.write(f"{family_lineage},{n_at_rank},{n_virus}\n")
            # currently not writing the actual COUNT of each unique one. not sure if we want/need that

# Command line interface
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate lineage information to final annotations.")
    parser.add_argument("all_idents_csv", help="CSV file with all identifiers for annotation")
    parser.add_argument("--dna-lineages", help="CSV file with dna-based lineage information")
    parser.add_argument("--protein-lineages", help="CSV file with protein-based lineage information")
    parser.add_argument("--dna-cluster-lineages", help="CSV file with dna-based annotated clusters")
    parser.add_argument("--protein-cluster-lineages", help="CSV file with protein-based annotated clusters")
    parser.add_argument("-o", "--output-csv", help="Output CSV file to save the results")
    parser.add_argument('-r', '--rank', default='species', help= 'rank to report counts')
    args = parser.parse_args()
    main(args)