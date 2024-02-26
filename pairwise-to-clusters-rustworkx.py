import argparse
import rustworkx as rx
from tqdm import tqdm
import csv
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

def plot_histogram(cluster_sizes, filename):
        fig, ax = plt.subplots()
        sns.set_theme(style="whitegrid")
        # Find the range of the data
        data_range = [np.min(cluster_sizes), np.max(cluster_sizes)]
        # Create the bins - one for each possible value in the range of data
        bins = np.arange(data_range[0], data_range[1] + 2) - 0.5
        sns.histplot(cluster_sizes, bins=bins, color='navy', kde=False, ax=ax)
        ax.set(xlabel='Cluster Size', ylabel='Frequency (Log Scale)', 
            title='Histogram of Cluster Sizes')

        # Use a logarithmic scale on y-axis to handle outliers
        ax.set_yscale('log')

        # Let matplotlib handle the x-axis ticks automatically
        plt.xticks(rotation=90)

        # Remove top and right borders
        sns.despine()

        # Save the plot as a png file
        fig.savefig(filename)
        plt.close()

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

            # Record only if similarity is above the threshold
            if similarity > threshold:

                for name in [query_name, match_name]:
                    if name not in name_to_id:
                        name_to_id[name] = current_id
                        nodes_to_add.add(current_id)
                        current_id += 1

                query_id = name_to_id[query_name]
                match_id = name_to_id[match_name]

                edges.append((query_id, match_id, similarity))
                if len(edges) >= batch_len:
                    graph.add_nodes_from(list(nodes_to_add))
                    nodes_to_add.clear()
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
    cluster_sizes = []
    with open(output_csv, 'w') as clusters_file:
        for component in connected_components:
            named_component = [id_to_name[node] for node in component]
            cluster_sizes.append(len(component))
            clusters_file.write(','.join(named_component) + '\n')

    print(f"Clusters written to {output_csv}")

    return connected_components, cluster_sizes

def main(args):
    print("Constructing graph...")
    graph, name_to_id = construct_graph(args.pairwise_csv, args.threshold, args.similarity_column)
    
    print("Clustering...")
    cc, cluster_sizes = cluster_graph(graph, args.output_csv, name_to_id)
    
    plot_histogram(cluster_sizes, args.histogram)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description="Constructs a graph and clusters sequences based on given parameters.")
    p.add_argument('pairwise_csv', help="Output from branchwater pairwise")
    p.add_argument('-t', '--threshold', required=True, type=float, help="distance threshold")
    p.add_argument('-d', '--similarity-column', required=True, help="similarity measure to use for clustering", choices=["ANI", "containment"])
    p.add_argument('-o', '--output_csv', required=True, help="CSV file to write clusters to")
    p.add_argument('-g', '--histogram', help="cluster size histogram")
    args = p.parse_args()
    main(args)