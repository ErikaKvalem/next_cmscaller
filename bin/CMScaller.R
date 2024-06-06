#!/usr/bin/env Rscript
'CMScaller.R
Usage:
  CMScaller.R --counts=<counts>

Mandatory arguments:
  --counts=<counts>                   count mat 

' -> doc

# Libraries
library(docopt)
library(biomaRt)
library(Biobase)
library(BiocStyle)
library(CMScaller)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(dplyr)
knitr::opts_chunk$set(fig.width=6, fig.height=3, 
                      dev.args=list(pointsize=8), dpi=150,
                      collapse=TRUE, message=TRUE, echo=TRUE, warnings=FALSE)
options(scipen=-1, digits=2)

# Input from command line
arguments <- docopt(doc, version = "0.1")
counts <- arguments$counts
counts <- read.table(counts, sep = ",", header=1)


#Test
#counts <- read.table("/data/scratch/kvalem/projects/2022/OLD/next_cmscaller/work/8f/fd9f8e6a811f97f565ce6f7e629079/counts.tsv", sep = ",", header=1)

colnames(counts)[1] <- "gene_name"
counts$ensembl_id = mapIds(org.Hs.eg.db,
                           keys=counts$gene_name, 
                           column="ENSEMBL",
                           keytype="SYMBOL",
                           multiVals="first")


rownames(counts) <- counts$gene_name
counts<-counts[,-grep("ensembl_id",colnames(counts))]
counts_mat <- counts[,-1]


# prediction and gene set analysis
par(mfrow=c(1,2))
res <- CMScaller(emat=counts_mat, RNAseq=TRUE, FDR=0.05, rowNames = "symbol")
table(pred=res$prediction)
head(res, n=3)

write.csv("CMScaller_results.csv")
