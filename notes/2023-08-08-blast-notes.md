new blast classifications are up in the google drive folder!

differences:
- used 2023-08-07 spillover file
- ICTV VMR database included nearly all entries for vmr 38. There are still two missing (Salmonella phage Fels2, Caenorhabditis elegans Cer13 virus). These are prophage that need to be extracted from their host genomes.

There are only 19 spillover accessions which had no hits via either blastn or blastx. They do not seem to correspond to those missing 2 reference VMRs. Instead, it looks like most of them were originally labeled as "Simian immunodeficiency virus", which we do have a VMR genome for. I ended up digging into these a bit, and found a useful illustrative example for thresholding.

When I looked in the new spillover file, it turns out there are 17,201 entries for this virus! So I looked at our blast results, and we seem to have matched most of them to four species in the same genus (Lentivirus):

```
1 Bovine immunodeficiency virus
848 Human immunodeficiency virus 1
 15626 Human immunodeficiency virus 2
709 Simian immunodeficiency virus
```
The bovine match was 97% pident but only for 38 base pairs length. The spillover accession came from a colobus monkey (https://www.ncbi.nlm.nih.gov/nuccore/AM937062.1/). The bovine VMR genome was the only blastn match, but there were a number of blastx matches to the simian immunodeficiency virus VMR genome, illustrating the need for blast thresholding to allow preferring the blastx match when the blastn isn't very good. The blastx matches were 57% match or lower, but over a much longer aligned region.

Calculating/thresholding on the aligned fraction (aln length/total length) for each match should prove useful here. I'll work on it this week, along with thresholds for blastx and thinking about when to pull matches up to e.g. genus level.