#!/usr/bin/env Rscript
'CMScaller.R
Usage:
  CMScaller.R --count_mat=<count_mat>

Mandatory arguments:
  --count_mat=<count_mat>                   count mat 

' -> doc

# Libraries
library(docopt)
library(Biobase)
library(BiocStyle)
library(CMScaller)
library("AnnotationDbi")
library("org.Hs.eg.db")
knitr::opts_chunk$set(fig.width=6, fig.height=3, 
                      dev.args=list(pointsize=8), dpi=150,
                      collapse=TRUE, message=TRUE, echo=TRUE, warnings=FALSE)
options(scipen=-1, digits=2)


# Input from command line
arguments <- docopt(doc, version = "0.1")
counts <- read.csv(arguments$count_mat, sep = "\t")
rownames(counts) <- counts$X
counts <- counts[,-1]

# prediction and gene set analysis
par(mfrow=c(1,2))
res <- CMScaller(emat=counts, RNAseq=TRUE, FDR=0.05, rowNames = "symbol")
table(pred=res$prediction)
head(res, n=3)

write.csv(res, "CMScaller_results.csv")