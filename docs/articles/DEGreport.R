## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE,
               fig.width=9,fig.height=5,
               message=FALSE)

## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle::latex()

## ----package-load,message=FALSE--------------------------------------------
library(DEGreport)
data(humanGender)

## ----chunk-1---------------------------------------------------------------
library(DESeq2)
idx <- c(1:10, 75:85)
dds <- DESeqDataSetFromMatrix(assays(humanGender)[[1]][1:1000, idx],
                              colData(humanGender)[idx,], design=~group)
dds <- DESeq(dds)
res <- results(dds)

## ----chunk-2---------------------------------------------------------------
counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))

## ----chunk-size-factor-----------------------------------------------------
degCheckFactors(counts[, 1:6])

## ----chunk-qc--------------------------------------------------------------
degQC(res[["pvalue"]], counts, design[["group"]])

## ----chunk-covariates------------------------------------------------------
resCov <- degCovariates(log2(counts(dds)+0.5),
                        colData(dds))
resCov$plot

