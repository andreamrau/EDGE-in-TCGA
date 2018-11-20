#---------------------------------------------------------------------------------------------------
# Load packages
library(prodlim)
library(survival)
library(survcomp)
library(BiocGenerics)
library(Biobase)
library(mclust)
library(limma)
library(AIMS)
library(biomaRt)
library(dplyr)
library(lazyeval)

## Data output from https://github.com/andreamrau/EDGE-in-TCGA/1_download_TCGA.R
load("BRCA/TCGA_formattedData.RData")

#----------------------------------------------------------------------------------------------------
# Subset data

## Only tumor samples, only data for which we have clinical information
tumor_index <- which(substr(colnames(rnaseq), 14, 15)  < 10)
rnaseq_tumor <- rnaseq[,tumor_index]
clinical_index <- which(substr(colnames(rnaseq_tumor), 1, 12) %in% rownames(clinical))
rnaseq_tumor_clinical <- rnaseq_tumor[,clinical_index]
dim(rnaseq_tumor_clinical) == dim(rnaseq_tumor)
# First need to identify the Entrez IDs for each gene (leave in duplicates for AIMS)
date <- "aug2017" ## Ensembl 90
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                   host=paste(date,"archive.ensembl.org",sep="."),
                   dataset="hsapiens_gene_ensembl")
entrez_ids <- getBM(attributes = c('entrezgene', 'hgnc_symbol'),
                    filters = 'hgnc_symbol',
                    values = rownames(rnaseq_tumor),
                    mart = ensembl)

rnaseq_df <- data.frame(hgnc_symbol = rownames(rnaseq_tumor), rnaseq_tumor, check.names=FALSE)
rnaseq_df <- rnaseq_df %>% left_join(., entrez_ids, by = "hgnc_symbol")
rnaseq_final <- as.matrix(dplyr::select(rnaseq_df, -entrezgene, -hgnc_symbol))
rownames(rnaseq_final) <- rnaseq_df$hgnc_symbol
entrez_final <- as.character(rnaseq_df$entrezgene)

#----------------------------------------------------------------------------------------------------
# AIMS subtypes

aims <- applyAIMS(rnaseq_final, entrez_final)
aims_subtypes <- data.frame(ID=substr(rownames(aims$cl), 1, 12), subtype=aims$cl[,1])
table(aims_subtypes$subtype)
## Remove metastatic samples
remove_metastatic <- grep("06A", rownames(aims_subtypes))
aims_subtypes <- aims_subtypes[-remove_metastatic,]
write.table(aims_subtypes, file="aims_subtypes.txt", col.names=FALSE, row.names=FALSE, sep="\t")
