#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate base 

awk 'NR!=1{print $1}' ../data/nextstrain_ncov_open_reference_metadata.tsv > ../results/Genbank_lineage_representative.list 

xzcat /Volumes/SSD_480G/Downloads/aligned_20220317.fasta.xz | seqtk subseq - ../results/Genbank_lineage_representative.list > ../results/Genbank_lineage_representative.fasta

conda deactivate
