import argparse
import rustworkx as rx
from tqdm import tqdm
import csv

def ani_from_containment(containment, ksize):
    if containment == 0.0:
        return 1.0
    elif containment == 1.0:
        return 0.0
    else:
        return 1.0 - (1.0 - containment ** (1.0 / ksize))

def construct_graph(pairwise_file, threshold, similarity_column="ANI", ksize=31, batch_len=10000):
    graph = rx.PyGraph()
    name_to_id = {}
    current_id = 0
    nodes_to_add = set()

    with open(pairwise_file, newline='') as f:
        reader = csv.DictReader(f, delimiter=',')
        # Check if 'ANI' is in the header or if we need to calculate it
        use_ani_conversion = similarity_column == "ANI" and similarity_column not in reader.fieldnames
        
        edges = []
        for row in tqdm(reader, desc="Processing rows"):
            query_name = row["query_name"]
            match_name = row["match_name"]
            if query_name == match_name:
                continue

            # Determine the similarity based on available data
            if use_ani_conversion:
                # Assume containment is available and convert it
                containment = float(row["containment"])
                similarity = ani_from_containment(containment, ksize) * 100.0
            else:
                similarity = float(row[similarity_column]) * 100.0

            # Continue only if similarity is above the threshold
            if similarity < threshold:
                continue

            for name in [query_name, match_name]:
                if name not in name_to_id:
                    name_to_id[name] = current_id
                    nodes_to_add.add(current_id)
                    current_id += 1

            if len(nodes_to_add) > batch_len:
                for node_id in nodes_to_add:
                    graph.add_node(node_id)
                nodes_to_add.clear()

            query_id = name_to_id[query_name]
            match_id = name_to_id[match_name]

            edges.append((query_id, match_id, similarity))
            if len(edges) >= batch_len:
                graph.add_edges_from(edges)
                edges.clear()

    if nodes_to_add:
        graph.add_nodes_from(list(nodes_to_add))
    if edges:
        graph.add_edges_from(edges)
            
    return graph, name_to_id
 
def cluster_graph(graph, output_csv, name_to_id):
    connected_components = rx.connected_components(graph)
    print(f"Number of clusters: {len(connected_components)}")
    
    id_to_name = {v: k for k, v in name_to_id.items()}
    with open(output_csv, 'w') as clusters_file:
        for component in connected_components:
            named_component = [id_to_name[node] for node in component]
            clusters_file.write(','.join(named_component) + '\n')

    print(f"Clusters written to {output_csv}")

def main(args):
    print("Constructing graph...")
    graph, name_to_id = construct_graph(args.pairwise_csv, args.threshold, args.similarity_column)
    
    print("Clustering...")
    cluster_graph(graph, args.output_csv, name_to_id)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="Constructs a graph and clusters sequences based on given parameters.")
    p.add_argument('pairwise_csv', help="Output from branchwater pairwise")
    p.add_argument('-t', '--threshold', required=True, type=float, help="distance threshold")
    p.add_argument('-d', '--similarity-column', required=True, help="similarity measure to use for clustering", choices=["ANI", "containment"])
    p.add_argument('-o', '--output_csv', required=True, help="CSV file to write clusters to")
    args = p.parse_args()
    main(args)