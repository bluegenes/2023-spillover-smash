import csv
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import argparse

def read_coords_file(coords_file):
    reference_mappings = {}
    with open(coords_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            genome, contig, start, end = row
            start, end = int(start), int(end)
            length = end - start + 1
            contig_name_with_coords = f"{contig}:{start}-{end}"
            reference_mappings[contig_name_with_coords] = (genome, contig, length)
    return reference_mappings

def closest_references_for_query(query_idx, distance_matrix, ref_indices, n):
    distances = [(ref_idx, distance_matrix[query_idx][ref_idx]) for ref_idx in ref_indices]
    distances.sort(key=lambda x: x[1])
    query_res = []
    for i, (r_idx, dist) in enumerate(distances):
        ref_name = distance_matrix.names[r_idx]
        if i == n:
            return query_res 
        query_res.append((ref_name, dist))
    return query_res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find closest reference in an MSA')
    parser.add_argument('-m', '--msa-fasta', type=str, help='MSA file in FASTA format')
    parser.add_argument('-c', '--coords-csv', type=str, help='Coordinates CSV file with reference genome mapping')
    parser.add_argument('-o', '--output-csv', type=str, help='Output CSV file')
    parser.add_argument("--top-n", type=int, default=2, help="Number of top references to output for each query.")

    args = parser.parse_args()

    msa = AlignIO.read(args.msa_fasta, 'fasta')
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(msa)

 # Index references and queries
    ref_mapping = read_coords_file(args.coords_csv)
    ref_indices = [distance_matrix.names.index(ref) for ref in ref_mapping.keys() if ref in distance_matrix.names]

    with open(args.output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Query', 'Reference Contig', 'Similarity', 'Alignment Length', 'GenBank Assembly ID'])
        
        for idx, name in enumerate(distance_matrix.names):
            if idx not in ref_indices:
                closest_refinfo = closest_references_for_query(idx, distance_matrix, ref_indices, args.top_n)
                for ref_contig, dist in closest_refinfo:
                    ref_genome, contig, length = ref_mapping[ref_contig]
                    writer.writerow([name,
                                     contig.rsplit(':')[0],
                                     1-dist,
                                     length,
                                     ref_genome])



# def closest_reference_for_query(query_idx, distance_matrix, ref_indices):
#     min_distance = float('inf')
#     closest_ref_idx = None

#     for ref_idx in ref_indices:
#         distance = distance_matrix[query_idx][ref_idx]
#         if distance < min_distance:
#             min_distance = distance
#             closest_ref_idx = ref_idx
    
#      # Extract coordinates from the reference name
#     closest_ref_name, coords = distance_matrix.names[closest_ref_idx].split(':')
#     start, end = map(int, coords.split('-'))
#     length = end - start + 1
    
#     return closest_ref_name, 1 - min_distance, length