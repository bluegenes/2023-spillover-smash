import csv
import networkx as nx
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree, fcluster
from scipy.spatial.distance import squareform
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import argparse

def build_newick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        node_name = leaf_names[node.id]
        return "%s:%.2f%s" % (node_name, parentdist - node.dist, newick)
    newick = (
        "):%.2f%s" % (parentdist - node.dist, newick)
        if len(newick) > 0
        else ");"
    )
    newick = build_newick(node.get_left(), newick, node.dist, leaf_names)
    newick = build_newick(node.get_right(), f",{newick}", node.dist, leaf_names)
    newick = f"({newick}"
    return newick


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


# convert graph to newick format
def graph_to_newick(filename, metric):
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

    # convert linkage matrix to tree
    print("Converting linkage matrix to tree...")
    tree = to_tree(link_matrix, rd=False)

    print("Building newick tree...")
    leaf_names = list(graph.nodes())
    newick = build_newick(tree, "", tree.dist, leaf_names)
    return newick

def main():

    parser = argparse.ArgumentParser(description='Convert pairwise file to newick format')
    parser.add_argument('--results', type=str, required=True, help='results csv')
    parser.add_argument('--metric', type=str, required=False, default = 'jaccard', help='metric')
    parser.add_argument('--output', type=str, required=True, help='Output file')
    args = parser.parse_args()


    pairwise_file = args.results
    metric = args.metric
    print(f"using metric '{metric}'")

    newick = graph_to_newick(pairwise_file, metric) # id_to_name


    with open(args.output, 'w') as OUT:
        OUT.write(newick)

if __name__ == "__main__":
    main()
