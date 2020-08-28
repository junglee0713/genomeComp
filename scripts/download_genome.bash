#!/bin/bash
taxa_assembly=$1
whole_dir=$2
rna_dir=$3
while IFS=$'\t' read taxa assembly; do
    taxa="${taxa// /_}"
    whole="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/${taxa}/latest_assembly_versions/${assembly}/${assembly}_genomic.fna.gz"
    wget ${whole} -P "${whole_dir}"
    rna="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/${taxa}/latest_assembly_versions/${assembly}/${assembly}_rna_from_genomic.fna.gz"
    wget ${rna} -P "${rna_dir}"
done < ${taxa_assembly}
