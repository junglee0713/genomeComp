import configparser
import yaml
import re

whole_genome_dir = config['all']['root']  + '/whole'
rna_dir = config['all']['root'] + '/rna'
ssu_dir = config['all']['root'] + '/16S'
ssu_blastDB = config['all']['root'] + '/16S_blastDB'
ANIb_out_dir = config['all']['root'] + '/ANIb_out'
taxa_assembly_fp = config['all']['taxa_assembly']

workdir: config['all']['root']

#######
#  create the list of assembly names
#######

with open(taxa_assembly_fp) as f:
    lines = f.readlines()
assembly_names = []
for line in lines:
    assembly_names.append(line.split('\t')[1].rstrip())

rule all:
    input:
        ANI_result = ANIb_out_dir + '/ANIb_percentage_identity.tab', 
        pctid_result = config['all']['root'] + '/blast_hits.tsv'

rule run_blast:
    input:
        sentinel = expand(ssu_blastDB + '/16S.fasta.{suffix}', suffix = ['nhr', 'nin', 'nsq']),
        db = ssu_blastDB + '/16S.fasta',
        query = ssu_blastDB + '/16S.fasta'
    output:
        config['all']['root'] + '/blast_hits.tsv'
    shell:
        """
            blastn -evalue 1e-5 -outfmt 7 \
            -db {input.db} \
            -query {input.query} \
            -num_threads 4 \
            -out {output}
        """

rule make_blastDB:
    input: 
        ssu_blastDB + '/16S.fasta'
    output:
        expand(ssu_blastDB + '/16S.fasta.{suffix}', suffix = ['nhr', 'nin', 'nsq'])
    shell:
        """
            makeblastdb -dbtype nucl -in {input}
        """

rule combine_16S:
    input:
        expand(ssu_dir + '/{assembly}_16S.fna', assembly = assembly_names)
    output:
        ssu_blastDB + '/16S.fasta'
    shell:
        """
            cat {input} > {output}
        """

rule extract_16S:
    input:
        rna_dir + '/{assembly}_rna_from_genomic.fna'
    output:
        ssu_dir + '/{assembly}_16S.fna'
    params:
        ssu_dir
    shell:
        """
            mkdir -p {params}
            okfasta searchdesc \
                --input {input} \
                --output {output} \
                "product=16S ribosomal RNA"
        """

rule ANIb:
    input:
        unzipped = expand(whole_genome_dir + '/{assembly}_genomic.fna', assembly = assembly_names)
    output:
        ANIb_out_dir + '/ANIb_percentage_identity.tab'
    params:
        indir = whole_genome_dir,
        outdir = ANIb_out_dir
    shell:
        """
            average_nucleotide_identity.py \
                -i {params.indir} \
                -o {params.outdir} \
                -m ANIb -g --force
        """

rule unzip_rna:
    input:
        rna_dir + '/{assembly}_rna_from_genomic.fna.gz'
    output:
        rna_dir + '/{assembly}_rna_from_genomic.fna'
    shell:
        """
            gunzip {input}
        """

rule unzip_whole:
    input:
        whole_genome_dir + '/{assembly}_genomic.fna.gz'
    output:
        whole_genome_dir + '/{assembly}_genomic.fna'
    shell:
        """
            gunzip {input}
        """

rule download_genome:
    input:
        taxa_assembly_fp
    output:
        whole = expand(whole_genome_dir + '/{assembly}_genomic.fna.gz', assembly = assembly_names),
        rna = expand(rna_dir + '/{assembly}_rna_from_genomic.fna.gz', assembly = assembly_names)
    params:
        script = config['download_genome']['script'],
        whole = whole_genome_dir,
        rna = rna_dir
    shell:
        """
            mkdir -p {params.whole} {params.rna}
            {params.script} {input} {params.whole} {params.rna}
        """

onsuccess:
    print('Workflow finished, no error')
    shell('mail -s "Workflow finished successfully" ' + config['all']['admin_email'] + ' < {log}')

onerror:
    print('An error occurred')
    shell('mail -s "An error occurred" ' + config['all']['admin_email'] + ' < {log}')

