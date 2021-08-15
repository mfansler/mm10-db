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
    input=list(),
    output=list(gtf="txdb.mm10.ensGene.txcutr.w{width}.gtf",
                fa="txdb.mm10.ensGene.txcutr.w{width}.fa"),
    wildcards=list(width="500"),
    threads=1
  )
}

################################################################################
## Libraries and Parameters
################################################################################

library(txcutr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

mm10 <- BSgenome.Mmusculus.UCSC.mm10
txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

maxTxLength <- as.integer(snakemake@wildcards$width)

## set cores
BiocParallel::register(BiocParallel::MulticoreParam(snakemake@threads))

################################################################################
## Truncate and Export
################################################################################

txdb_result <- truncateTxome(txdb, maxTxLength)

exportGTF(txdb_result, snakemake@output$gtf)

exportFASTA(txdb_result, mm10, snakemake@output$fa)
