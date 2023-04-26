# 2023-spillover-smash

This is a workflow to reclassify the SpillOver sequences using all ICTV genomes. It relies on snakemake for automation and conda for software management. 

## Build ICTV Database

### Find Corresponding 'Assembly' datasets

The ICTV viruses are provided via an excel file with one row per reference virus, and contains taxonomic information as well as additional reference information. This file provides NCBI accession numbers in GenBank and RefSeq that correspond to the sequences. For Segmeneted Viruses, an accession number is provided per genome segment, sometimes labeled, e.g. `{name}: {accession}`. For classification of the entire virus, it would be ideal instead to have a combined accession that refers to all genomic sequence of each virus. NCBI provides this in their curated 'Assembly' datsets resource, and we can find the corresponding identifier by querying via NCBI Entrez or using the `ncbi-datasets` tool. 

The file `find-assembly-accessions` uses the `Biopython Entrez` library and produces an updated `VMR` file with the corresponding GenBank Assembly accession. There are a few exceptions where either no corresponding Assembly exists, or it was not able to be matched via entrez.

Exceptions:
1. No Accession was provided (`no_accession`; 232)
    - No GenBank accession was provided, so we don't have a sequence record to link

2. Viral genomes within other genomes (`parentheses`; 2)
    - The accessions provided here correspond to larger genomes, where the viral genome is a small component, with base pair range provided in parantheses. Since we don't want to assign the entire host genome to a viral taxonomy, we should manually download the base pair range and use it as the sequence. Alternatively, in some cases there may be a better representative for this species available.

3. No Assembly Accession was found (`no_assembly`; 154): 
    - Partial genomes (54)
    - ?? (100)

4. Multiple Assembly Accessions found (`multiple_acc`; 13)
    - The segments provided for these viruses link to different Assembly accessions. This should not happen. In the one case I spot-checked, this was due to an assembly being redacted. Needs manual investigation to select the best accession.

5. Error during information retrieval (`retrieval`; 0)
    - Primarily a check for any errors occurring during entrez database queries.


### Download Viral Assemblies and build a `sourmash` database

With assembly accessions in hand, we now download all viral genomes from NCBI and index them as a sourmash database. We use the file generated above as input into the `dl-sketch.ICTV.smk` snakefile. 

This is currently written to download and sketch all available assembly accessions, but we may want to build a smaller database containing only the `exemplar` genomes, leaving out the `additional` isolates. In the `VMR_21-221122_MSL37.csv` file, there are 10435 Exemplar and 1626 'Additional'.

Issues:
- Ran into an exception during download: GCA_004789135.1 is supressed. For now, manually removed from accession file. In the future, need programmatic query for determining whether record has been supressed -- should be feasible to add to the first step, where we use biopython entrez to link assembly identifiers.
- 6 records have duplicate entries in the csv file:
  > - GCF_001041915.1 Mycobacterium phage Fionnbharth
  > - GCF_001745335.1 Shigella phage SHBML-50-1
  > - GCF_002625825.1 Rhodococcus phage Hiro
  > - GCF_003308095.1 Rhodococcus phage Takoda
  > -GCF_004138835.1 Salmonella phage 3-29
  > - GCF_004138895.1 Pseudomonas phage vB_PaeM_SCUT-S1
  
  > Since the duplicated entries in the `fromfile.csv` file are identical, they were removed via `awk '!x[$0]++' $filename > $filename-new` for now. This should be automated in the future.


## Build SpillOver Database

For the SpillOver sequences, we can download the provided accessions directly. We use the `2023-03-27_spillover_accession-numers.csv` as input into `dl-sketch.spillover.smk` snakefile, selecting accessions from the `AccessionNumber` column. 

Exceptions:
- Some entries (384 of 34220 total) contain `Unknown` accessions, which means we cannot download their sequences for reclassification.
- 269 were duplicated entries. As above, these were removed via `awk '!x[$0]++' $filename > $filename-new` for now, checked, and then moved back to the original filename. This should be automated in the future.

## Classify SpillOver Sequences
