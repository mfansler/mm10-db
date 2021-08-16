#!/usr/bin/env Rscript

################################################################################
## Mock `snakemake` preamble
################################################################################

if (interactive()) {
  library(methods)
  Snakemake <- setClass(
    "Snakemake", 
    slots=c(
      input='list', 
      output='list',
      wildcards='list',
      threads='numeric'
    )
  )
  snakemake <- Snakemake(
      input=list(gtf="gencode.vM21.annotation.mRNA_ends_found.gtf.gz"),
      output=list(gtf="/fscratch/fanslerm/gencode.vM21.annotation.mRNA_ends_found.txcutr.w500.gtf",
                  fa="/fscratch/fanslerm/gencode.vM21.annotation.mRNA_ends_found.txcutr.w500.fa"),
      wildcards=list(width="500"),
      threads=1
  )
}

################################################################################
## Libraries and Parameters
################################################################################

library(txcutr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicFeatures)

mm10 <- BSgenome.Mmusculus.UCSC.mm10

maxTxLength <- as.integer(snakemake@wildcards$width)

## set cores
BiocParallel::register(BiocParallel::MulticoreParam(snakemake@threads))

################################################################################
## Load Data, Truncate, and Export
################################################################################

txdb <- makeTxDbFromGFF(file=snakemake@input$gtf,
                        organism="Mus musculus",
                        taxonomyId=10090,
                        chrominfo=seqinfo(mm10))

txdb_result <- truncateTxome(txdb, maxTxLength)

exportGTF(txdb_result, snakemake@output$gtf)

exportFASTA(txdb_result, mm10, snakemake@output$fa)
