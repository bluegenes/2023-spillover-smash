import sourmash
import pandas as pd
import numpy as np 

from sourmash.tax.tax_utils import ICTVRankLineageInfo


# make ICTVLineageInfo object for each lineage
def make_ICTVLineageInfo(row, column='blast_lineage'):
    lin = row[column]
    if not isinstance(lin, str): # if nan
        return ICTVRankLineageInfo()
    lin_str = lin.split(';', 1)[1] # drop 'Viruses;'
    lineage =  ICTVRankLineageInfo(lineage_str = lin_str)
    return lineage

# extract lineage to columns
def lin_to_cols(row, column='blast_linfo', ranks_to_keep = ['name','species', 'genus', 'family']):
    prefix = 'blast_vmr_'
    if 'gather' in column:
        prefix = 'gather_vmr_'
    lin = row[column]
    lin_dict = {prefix + a.rank: a.name for a in lin.lineage}
    for k, v in lin_dict.items():
        if k.split(prefix)[1] in ranks_to_keep:
            row[k] = v
    return row


def main(args):
    # read in metadata info
    sDF = pd.read_csv(args.info)
    blastDF = pd.read_csv(args.blast, sep = '\t')

    # make sourmash taxonomy out of lineage
    blastDF['blast_linfo'] = cDF.apply(make_ICTVLineageInfo, axis=1, column='taxonomy')
    blastDF = blastDF.apply(lin_to_cols, axis=1, column='blast_linfo')

    merged = sDF.merge(blastDF, right_on='qseqid', left_on='AccessionNumber', how='left') # 'left' to keep all spillover accessions
    merged_write = merged.drop(columns=['blast_linfo'])
    merged_write.to_csv(args.output, sep = '\t', index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    parser = argparse.ArgumentParser(description="Assess classification")
    parser.add_argument("--info", type=str, help="Metadata for original accessions", default = '../inputs/2023-03-27_spillover_accession-numers.csv')
    parser.add_argument("--blast", type=str, help="Blast classifications")
    parser.add_argument("--output", type=str, help="output combined file")

    # Parse command line arguments
    args = parser.parse_args()

    # Download sequences for the specified accession
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
