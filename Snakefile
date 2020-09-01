#!/usr/bin/env snakemake

configfile: "config.yaml"

import os

# make sure the tmp directory exists
os.makedirs(config["tmp_dir"], exist_ok=True)


rule all:
    input: ""

rule blacklist:
    output:
        "mm10.blacklist.bed.gz"
    params:
        url="http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz"
    shell:
        """
        wget -O {output} {params.url}
        """

rule polyAsite:
    output:
        bed="polyASite.clusters.mm10.bed.gz",
        bed_tissue="polyASite.clusters.tissues.mm10.bed.gz",
        idx="polyASite.clusters.mm10.bed.gz.tbi",
        idx_tissue="polyASite.clusters.tissues.mm10.bed.gz.tbi"
    params:
        url="https://polyasite.unibas.ch/clusters/Mus_musculus/4.0/clusters.bed",
        url_tissue="https://polyasite.unibas.ch/clusters/Mus_musculus/4.0/clusters_withTissueInfo.bed"
    version: 1.0
    shell:
        """
        wget -O- {params.url} | sort -k1,1 -k2,2n | bgzip > {output.bed}
        tabix {output.bed}
        wget -O- {params.url_tissue} | sort -k1,1 -k2,2n | bgzip > {output.bed_tissue}
        tabix {output.bed_tissue}
        """

rule gencode:
    output:
        gff="gencode.vM{version}.annotation.gff3.gz",
        gtf="gencode.vM{version}.annotation.gtf.gz",
        fa="gencode.vM{version}.transcripts.fa.gz"
    params:
        url_base = lambda wcs: "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M%s/" % wcs.version
    wildcard_constraints:
        version = '\d+'
    shell:
        """
        wget {params.url_base}{output.gff}
        wget {params.url_base}{output.gtf}
        wget {params.url_base}{output.fa}
        """

rule clean_gencode_fasta:
    input:
        "gencode.vM{version}.{collection}.fa.gz"
    output:
        "gencode.vM{version}.{collection}.clean.fa.gz"
    wildcard_constraints:
        version = '\d+'
    threads: 4
    shell:
        """
        gzip -cd {input} | sed 's/|.*//' | pigz -p {threads} > {output}
        """

rule ucsc_fasta:
    output:
        "ucsc.mm10.fa"
    shell:
        """
        wget -qO - 'http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz' |
        gunzip -c > {output}
        """

rule kallisto_index:
    input:
        "gencode.vM{version}.{collection}.clean.fa.gz"
    output:
        "gencode.vM{version}.{collection}.kdx"
    conda:
        "envs/kallisto.yaml"
    resources:
        mem=16
    shell:
        """
        kallisto index -i {output} {input}
        """

rule sort_gff:
    input:
        "gencode.vM{version}.annotation.gff3.gz"
    output:
        gff="gencode.vM{version}.annotation.sorted.gff3.gz",
        idx="gencode.vM{version}.annotation.sorted.gff3.gz.tbi"
    threads: 4
    shell:
        """
        gzip -cd {input} | awk '!( $0 ~ /^#/ )' | sort --parallel={threads} -S4G -k1,1 -k4,4n | bgzip -c > {output.gff}
        tabix {output.gff}
        """

rule refseq_to_symbol:
    output:
        "refseq2gene.tsv"
    shell:
        """
        mysql --user=genome -N --host=genome-mysql.cse.ucsc.edu -A -D mm10 -e 'select name,name2 from refGene' > {output}
        """

rule clean_gencode_gff_mRNA_ends:
    input:
        "gencode.vM{version}.annotation.sorted.gff3.gz"
    output:
        gff="gencode.vM{version}.annotation.mRNA_ends_found.gff3.gz",
        idx="gencode.vM{version}.annotation.mRNA_ends_found.gff3.gz.tbi"
    shell:
        """
        gzip -cd {input} |
        awk '$0 !~ /mRNA_end_NF/' |
        bgzip -c > {output.gff}
        tabix {output.gff}
        """

rule clean_gencode_gtf_mRNA_ends:
    input:
        "gencode.vM{version}.annotation.gtf.gz"
    output:
        gtf="gencode.vM{version}.annotation.mRNA_ends_found.gtf.gz"
    shell:
        """
        gzip -cd {input} |
        awk '$0 !~ /mRNA_end_NF/' |
        bgzip -c > {output.gtf}
        """

rule gencode_gtf2bed:
    input:
        "gencode.vM{version}.{annot}.gtf.gz"
    output:
        "gencode.vM{version}.{annot}.bed"
    conda:
        "envs/gtf2bed.yaml"
    shell:
        """
        gzip -cd {input} |
        grep -v '^#' |
        gtfToGenePred /dev/stdin /dev/stdout |
        genePredToBed stdin {output}
        """

rule gtf_tx2gene:
    input:
        "gencode.vM{version}.{annot}.gtf.gz"
    output:
        "gencode.vM{version}.{annot}.tx2gene.tsv"
    shell:
        """
        gzip -cd {input} | awk -f scripts/gtf_tx2gene.awk > {output}
        """

rule qapa_bed:
    output:
        "qapa_3utrs.gencode_VM22.mm10.bed"
    shell:
        """
        wget -qO - 'https://github.com/morrislab/qapa/releases/download/v1.3.0/qapa_3utrs.gencode_VM22.mm10.bed.gz' |
        gunzip -c > {output}
        """

rule qapa_fasta:
    input:
        fa="ucsc.mm10.fa",
        bed="qapa_3utrs.gencode_VM22.mm10.bed"
    output:
        "qapa_3utrs.gencode_VM22.mm10.fa"
    singularity:
        "docker://brianyee/qapa:1.3.0"
    shell:
        """
        qapa fasta -f {input.fa} {input.bed} {output}
        """

rule qapa_salmon_idx:
    input:
        "qapa_3utrs.gencode_VM22.mm10.fa"
    output:
        directory("qapa_3utrs.gencode_VM22.mm10.sidx")
    conda: "envs/salmon.yaml"
    resources:
        mem=8
    shell:
        """
        salmon index -t {input} -i {output}
        """

rule qapa_salmon_decoy_idx:
    input:
        qapa=rules.qapa_fasta.output,
        genome=rules.ucsc_fasta.output
    params:
        fa=config['tmp_dir'] + "/mm10_plus_qapa.fasta",
        decoy=config['tmp_dir'] + "/decoys.txt"
    output:
        directory("qapa_3utrs.gencode_VM22.mm10.decoy.sidx")
    threads: 12
    resources:
        mem=3
    conda: "envs/salmon.yaml"
    shell:
        """
        grep "^>" {input.genome} | cut -d ' ' -f 1 > {params.decoy}
        sed -i.bak -e 's/>//g' {params.decoy}
        cat {input.qapa} {input.genome} > {params.fa}
        salmon index -t {params.fa} -d {params.decoy} -p {threads} -i {output}
        rm {params.decoy}
        rm {params.fa}
        """

rule ensembl_ids:
    output:
        "ensembl_identifiers_v{version}.txt"
    shell:
        """
        mysql --user anonymous --host=martdb.ensembl.org --port=5316 -A ensembl_mart_{wildcards.version} \\
        -e "select stable_id_1023 as 'Gene stable ID', stable_id_1066 as 'Transcript stable ID', \\
        biotype_1020 as 'Gene type', biotype_1064 as 'Transcript type', \\
        display_label_1074 as 'Gene name' from mmusculus_gene_ensembl__transcript__main" \\
        > {output}
        """
