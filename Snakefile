#!/usr/bin/env snakemake

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
