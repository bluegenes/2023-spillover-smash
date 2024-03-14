import os, sys
import csv
import argparse

import sourmash
from sourmash.tax import tax_utils
from sourmash.logging import notify, error



def main(args):
    """
    Annotate cluster results with LCA taxonomic lineage.
    """

    # load taxonomic assignments
    tax_assign = tax_utils.MultiLineageDB.load_from_gather_with_lineages(args.taxonomy_csv, keep_full_identifiers=False,
                                                               keep_identifier_versions=False, force=args.force,
                                                               lins=args.lins, ictv=args.ictv)

    if not tax_assign:
        error(
            f'ERROR: No taxonomic assignments loaded from {",".join(args.taxonomy_csv)}. Exiting.'
        )
        sys.exit(-1)

    # get csv from args
    input_csvs = tax_utils.collect_gather_csvs(args.cluster_csv)

    # handle each gather csv separately
    for n, in_csv in enumerate(input_csvs):
        # Check for a column we can use to find lineage information:
        with sourmash.FileInputCSV(in_csv) as r:
            header = r.fieldnames
            # check for empty file
            if not header:
                notify(f"Cannot read from '{in_csv}'. Is file empty?")
                sys.exit(-1)

            # check cluster file
            if 'nodes' not in header:
                notify(f"no 'nodes' column. Is this {in_csv} cluster file?")
                sys.exit(-1)

            # make output file for this input
            out_base = os.path.basename(in_csv.rsplit(".csv")[0])
            this_outfile = os.path.join(args.out_dir, out_base, ".with-lineages.csv")
            out_header = header + ["lineage"]

            with sourmash.FileOutputCSV(this_outfile) as out_fp:
                w = csv.DictWriter(out_fp, out_header)
                w.writeheader()

                n = 0
                n_missed = 0
                rows_missed = 0
                for n, row in enumerate(r):
                    # find annotation for each node in the cluster, then take LCA.
                    cluster_annot = set()
                    nodes = row['nodes'] # should be the list of queries in this cluster
                    for node in nodes:
                        # find lineage and write annotated row
                        lineage=None
                        lin = tax_assign.get(node)
                        if lin:
                            if args.lins:
                                lineage = tax_utils.LINLineageInfo()
                            elif args.ictv:
                                lineage = tax_utils.ICTVRankLineageInfo()
                            else:
                                lineage = tax_utils.RankLineageInfo()
                            # add match lineage to cluster_annot
                            cluster_annot.add(lineage)
                        else:
                            n_missed += 1

                    # get LCA of the node taxonomic assignments
                    if len(cluster_annot) == 0:
                        rows_missed +=1
                        continue
                    elif len(cluster_annot) > 1:
                        lin_tree = sourmash.tax.LineageTree().add_lineages(cluster_annot)
                        lca_lin = lin_tree.find_lca()
                    else:
                        lca_lin = cluster_annot[0]

                    # writerow
                    row["lineage"] = lca_lin.display_lineage()

                    # write row to output
                    w.writerow(row)

                rows_annotated = (n + 1) - rows_missed
                if n_missed:
                    notify(f"Missed {n_missed} taxonomic assignments during annotation.")
                if not rows_annotated:
                    notify(
                        f"Could not annotate any rows from '{in_csv}'."
                    )
                    sys.exit(-1)
                else:
                    notify(
                        f"Annotated {rows_annotated} of {n+1} total rows from '{in_csv}'."
                    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate cluster results with LCA taxonomic lineage.")
    parser.add_argument("cluster_csv", help="Path to the cluster CSV file.")
    parser.add_argument("-t", "--taxonomy_csv", help="Path to the taxonomy CSV file.")
    parser.add_argument("-o", "--output_dir", help="Output directory for annotated CSV files.", default=".")
    parser.add_argument("-f", "--force", help="Force loading taxonomic assignments.", action="store_true")
    parser.add_argument("--lins", help="Use LIN taxonomic assignments.", action="store_true")
    parser.add_argument("--ictv", help="Use ICTV taxonomic assignments.", action="store_true")
    args = parser.parse_args()
    main(args)
