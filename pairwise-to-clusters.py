import csv
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree, fcluster
from scipy.spatial.distance import squareform
import pandas as pd
import argparse

def load_edges(filename, metric="ochiai"):
    print(f"Loading edges from {filename}...")
    edges = []
    with open(filename, newline='') as f:
        reader = csv.DictReader(f, delimiter=',')
        for n, row in enumerate(reader):
            if n % 10000 == 0:
                print(f'loading row {n}...')
            query_name = row["query_name"]
            match_name = row["match_name"]
            if query_name == match_name:
                continue

            samples_pair = sorted([query_name, match_name])
            distance = float(row[metric])
            edges.append((samples_pair[0], samples_pair[1], distance))

    print(f"Loaded {len(edges)} comparisons.")
    return edges

def graph_to_clusters(filename, metric, threshold):
    edges = load_edges(filename, metric)
    graph = nx.Graph()
    graph.add_weighted_edges_from(edges)

    # get unique nodes
    nodes = list(graph.nodes())

    # create a pandas DataFrame for a distance matrix
    df = pd.DataFrame(data=np.zeros((len(nodes), len(nodes))), index=nodes, columns=nodes)

    # fill DataFrame with distances where edges exist (distance = 1 - similarity)
    print("Calculating distances...")
    for node1, node2, similarity in edges:
        df.at[node1, node2] = 1 - similarity
        df.at[node2, node1] = 1 - similarity

    # create 1D condensed distance matrix
    print("Creating condensed distance matrix...")
    condensed_matrix = squareform(df.values, checks=False)

    # create linkage matrix
    print("Creating linkage matrix...")
    link_matrix = linkage(condensed_matrix, method='ward')

    # generate clusters
    print("Generating clusters...")
    clusters = fcluster(link_matrix, threshold, criterion='distance')
    return clusters, nodes

def main():

    parser = argparse.ArgumentParser(description='Convert pairwise file to clusters')
    parser.add_argument('--results', type=str, required=True, help='results csv')
    parser.add_argument('--measure', type=str, required=False, default='jaccard', help='measure')
    parser.add_argument('-t', '--threshold', type=float, required=True, help='Similarity threshold')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
    args = parser.parse_args()

    pairwise_file = args.results
    measure = args.measure
    threshold = args.threshold
    print(f"using measure '{measure}' with threshold {threshold}")

    clusters, nodes = graph_to_clusters(pairwise_file, measure, threshold)

    # Write clusters to output file
    with open(args.output, 'w') as OUT:
        for node, cluster_id in zip(nodes, clusters):
            OUT.write(f"{node},{cluster_id}\n")

if __name__ == "__main__":
    main()
