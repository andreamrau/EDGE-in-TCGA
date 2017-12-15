#-----------------------------------------------------------------------------------
# Read in arguments from shell, load packages, load downloaded data

## Read in arguments from shell, if wanted
args <- commandArgs(trailingOnly=TRUE)
cancer <- as.character(args[1])
subtype <- as.character(args[2])
# cancer <- "THYM"
# subtype <- "all"

number_sPCs <- 5 
rerun_everything <- FALSE

cat("***", cancer, ":", subtype, "***\n")

user <- system("echo $USER", intern=TRUE)
if(user == "pauer") {
  libpath <- "/raid-04/SPH/raua/R/x86_64-redhat-linux-gnu-library/3.2"
  basedir <- paste0("/raid-04/SPH/raua/TCGA_data/", cancer, "/")
}
if(user == "raua")  {
  libpath <- "/home/raua/Data/R/x86_64-redhat-linux-gnu-library/3.2"
  basedir <- paste0("/home/raua/Data/TCGA_data/", cancer, "/")
}
library(mixOmics, lib=libpath)
library(dplyr, lib=libpath)
library(tidyr, lib=libpath)

## Load TCGA data
load(paste0(basedir, "TCGA_formattedData.RData"))

if(rerun_everything) {
#-----------------------------------------------------------------------------------
## Identify BRCA subtypes, if desired
if(cancer == "BRCA" & subtype != "all") {
  ## HER2+, HER2-, Triple-negative status (0=luminal A, 1=TN)
  ## Use the her2 column (0=negative, 1=positive)
  ## Use the tn_c column (0=luminal A, 1=triple negative)
  BRCA_subtypes <- read.csv("annotation/TCGA_BC_Subtypes.csv", stringsAsFactors=FALSE)
  if(subtype == "HER2pos") 
    subtype_barcodes <- BRCA_subtypes$bcr_pt_barcode[which(BRCA_subtypes$her2 == 1)]
  if(subtype == "HER2neg")
    subtype_barcodes <- BRCA_subtypes$bcr_pt_barcode[which(BRCA_subtypes$her2 == 0)]
  if(subtype == "luminalA") 
    subtype_barcodes <- BRCA_subtypes$bcr_pt_barcode[which(BRCA_subtypes$tn_c == 0)]
  if(subtype == "tripleneg")
    subtype_barcodes <- BRCA_subtypes$bcr_pt_barcode[which(BRCA_subtypes$tn_c == 1)]
}

## Combine mirna_ga and mirna_hiseq, if they both exist
## If both exist, we prefer to use mirna_hiseq
ga_keep_index <- NA
if(is.null(mirna_ga)) mirna <- mirna_hiseq
if(is.null(mirna_hiseq)) mirna <- mirna_ga
if(!is.null(mirna_ga) & !is.null(mirna_hiseq)) {
  ga_keep_index <- which(!substr(colnames(mirna_ga), 1, 16) %in% 
    substr(colnames(mirna_hiseq), 1, 16))
  if(!length(ga_keep_index)) mirna <- mirna_hiseq
  if(length(ga_keep_index)) {
    mirna <- cbind(mirna_hiseq, mirna_ga[,ga_keep_index])
    names(ga_keep_index) <- colnames(mirna_ga)[ga_keep_index]
  }
}

## Remove all samples arising from healthy tissue
## The 4th entry of TCGA barcodes corresponds to tissue type
## 01-09=tumor, 10-19=normal, 20-29=control
## NB: mutation data is per patient (no need to subset here)
rnaseq <- rnaseq[,which(as.numeric(substr(colnames(rnaseq), 14, 15)) < 10)]
cna <- cna[,which(as.numeric(substr(colnames(cna), 14, 15)) < 10)]
methyl <- methyl[,which(as.numeric(substr(colnames(methyl), 14, 15)) < 10)]
mirna <- mirna[,which(as.numeric(substr(colnames(mirna), 14, 15)) < 10)]

## Remove all samples not identified as white
clinical <- clinical[which(clinical[,"race"] == "white"),]

## Identify set of barcodes common to all data types, and subset data 
if(subtype=="all") {
  common_barcodes <- Reduce(intersect, 
    list(clinical=rownames(clinical),
       rnaseq=substr(colnames(rnaseq), 1, 12),
       cna=substr(colnames(cna), 1, 12),
       mut=colnames(mut),
       mirna=substr(colnames(mirna), 1, 12),
       methyl=substr(colnames(methyl), 1, 12)))
}
if(subtype != "all") {
  common_barcodes <- Reduce(intersect, 
    list(clinical=rownames(clinical),
       rnaseq=substr(colnames(rnaseq), 1, 12),
       cna=substr(colnames(cna), 1, 12),
       mut=colnames(mut),
       mirna=substr(colnames(mirna), 1, 12),
       methyl=substr(colnames(methyl), 1, 12),
       subtype_barcodes=subtype_barcodes))
}

       
mut_subset <- mut[,match(common_barcodes, colnames(mut))]
clinical_subset <- clinical[match(common_barcodes, rownames(clinical)),]
methyl_subset <- methyl[, match(common_barcodes, substr(colnames(methyl), 1, 12))]
cna_subset <- cna[, match(common_barcodes, substr(colnames(cna), 1, 12))]
rnaseq_subset <- rnaseq[, match(common_barcodes, substr(colnames(rnaseq), 1, 12))]
mirna_subset <- mirna[, match(common_barcodes, substr(colnames(mirna), 1, 12))]

dim(mut_subset)
dim(clinical_subset)
dim(methyl_subset)
dim(cna_subset)
dim(rnaseq_subset)
dim(mirna_subset)

#-----------------------------------------------------------------------------------
## Remove genes, mirs, CNA with all zeros, and mut with <= 1 total samples with mutations

RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

rnaseq_subset <- rnaseq_subset[which(rowSums(rnaseq_subset) > 0),]
mirna_subset <- mirna_subset[which(rowSums(mirna_subset) > 0),]
cna_subset <- cna_subset[which(RowVar(cna_subset) > 0),]
mut_subset <- mut_subset[which(rowSums(mut_subset) > 1),]

#-----------------------------------------------------------------------------------
## Choose most variable methylation probes in gene promoter regions
## Load methylation probe annotation (IlmnID <-> Gene Name <-> Grouping)
## Only retain prboes that coresspond to within 1500bp of a TSS (TSS200, TSS1500)
## For each gene, choose the most variable probe as representative

methyl_annot <- read.csv(paste0("/home/raua/Data/TCGA_data/annotation/",
  "HumanMethylation450_15017482_v1-2.csv"), stringsAsFactors=FALSE) %>% 
  as.data.frame() %>%
  select(IlmnID, UCSC_RefGene_Name, UCSC_RefGene_Group) %>%
  separate_rows(., UCSC_RefGene_Name, UCSC_RefGene_Group, sep=";") %>%
  filter(UCSC_RefGene_Group %in% c("TSS1500", "TSS200"))

methyl_pergene_subset <- matrix(NA, nrow=nrow(rnaseq_subset), ncol=ncol(rnaseq_subset))
colnames(methyl_pergene_subset) <- colnames(methyl_subset)
rownames(methyl_pergene_subset) <- rownames(rnaseq_subset)
for(i in seq_len(nrow(rnaseq_subset))) {
  if(i/1000 == ceiling(i/1000)) cat("Gene",i,"\n")
  gene_name <- rownames(rnaseq_subset)[i]
  probe_find <- grep(gene_name, methyl_annot[,"UCSC_RefGene_Name"])
  if(!length(probe_find)) next;
  if(length(probe_find) == 1) {
    index <- which(rownames(methyl_subset) == methyl_annot[probe_find, "IlmnID"])
    methyl_pergene_subset[i,] <- 
      methyl_subset[index,]
    methyl_annot <- methyl_annot[-probe_find,]
  }
  if(length(probe_find) > 1) {
    tmp <- methyl_subset[which(rownames(methyl_subset) %in% methyl_annot[probe_find, "IlmnID"]),]
    if(is.null(nrow(tmp))) methyl_pergene_subset[i,] <- tmp
    if(!is.null(nrow(tmp))) {
      var_calc <- RowVar(tmp)      
      var_choose <- which.max(var_calc)
      if(!length(var_choose)) methyl_pergene_subset[i,] <- NA
      if(length(var_choose)) methyl_pergene_subset[i,] <- tmp[var_choose,]
    }
    methyl_annot <- methyl_annot[-probe_find,]
  }
}

## Transform methylation data to logit scale (no true 0's or 1's)
## If there are any, replace true 0's or 1's with closest values 
methyl_pergene_subset <- pmax(methyl_pergene_subset, 
  min(methyl_pergene_subset[methyl_pergene_subset>0], na.rm=TRUE))
methyl_pergene_subset <- pmin(methyl_pergene_subset, 
  max(methyl_pergene_subset[methyl_pergene_subset<1], na.rm=TRUE))
methyl_pergene_subset <- log(methyl_pergene_subset / (1 - methyl_pergene_subset))


#-----------------------------------------------------------------------------------
## Identify TFs

IPA_TFs <- read.table(paste0("/home/raua/Data/TCGA_data/annotation/",
  "2016-12-22_IPA_TFlist.txt"), stringsAsFactors=FALSE)$V1
TRRUST_TFs <- read.table(paste0("/home/raua/Data/TCGA_data/annotation/",
  "2017-03-14_TRRUST_TFlist.txt"), stringsAsFactors=FALSE)$V1
TF_full_list <- unique(c(IPA_TFs, TRRUST_TFs))
TF_list <- TF_full_list[which(TF_full_list %in% rownames(rnaseq_subset))]  # 877


#-----------------------------------------------------------------------------------
## Match up mutation and CNA data to genes in RNA-seq

index <- na.omit(match(rownames(mut_subset), rownames(rnaseq_subset)))
remove_index_mut <- which(is.na(match(rownames(mut_subset), rownames(rnaseq_subset))))
mut_pergene_subset <- matrix(NA, nrow=nrow(rnaseq_subset), ncol=ncol(mut_subset))
colnames(mut_pergene_subset) <- colnames(mut_subset)
rownames(mut_pergene_subset) <- rownames(rnaseq_subset)
mut_pergene_subset[index,] <- mut_subset[-remove_index_mut,]

index <- na.omit(match(rownames(cna_subset), rownames(rnaseq_subset)))
remove_index_cna <- which(is.na(match(rownames(cna_subset), rownames(rnaseq_subset))))
cna_pergene_subset <- matrix(NA, nrow=nrow(rnaseq_subset), ncol=ncol(cna_subset))
colnames(cna_pergene_subset) <- colnames(cna_subset)
rownames(cna_pergene_subset) <- rownames(rnaseq_subset)
cna_pergene_subset[index,] <- as.matrix(cna_subset[-remove_index_cna,])



#-----------------------------------------------------------------------------------
## Calculate mean values for rnaseq_subset, mut_subset, methyl_pergene_subset, 
## cna_pergene_subset, as well as sample size for each data type

mean_rnaseq_subset <- rowMeans(rnaseq_subset, na.rm=TRUE)
mean_mut_subset <- rowMeans(mut_pergene_subset, na.rm=TRUE)
mean_methyl_subset <- rowMeans(methyl_pergene_subset, na.rm=TRUE)
mean_cna_subset <- rowMeans(cna_pergene_subset, na.rm=TRUE)

nmol <- c(nrow(rnaseq_subset), nrow(mirna_subset), 
  sum(!is.na(rowSums(cna_pergene_subset))),
  sum(!is.na(rowSums(methyl_pergene_subset))),
  sum(!is.na(rowSums(mut_pergene_subset))), ncol(rnaseq_subset))
names(nmol) <- c("rnaseq", "mirna", "cna", "methyl", "mut", "samples")

}
if(!rerun_everything) {
  load(paste0(basedir, cancer, "_results.RData"))
} 

#-----------------------------------------------------------------------------------
## Perform sPCA on log-transformed TF data
## One analysis performed over all TFs, and removing each TF one by one (877 + 1)
## Choose first 5 uncorrelated sparse principal components 
## Special case for LGG (due to full deletions of one arm of chr1 and one arm of chr19,
##    causing high multicollinearity between CNA data and first sPC for TFs):
##    remove the first sPC from consideration as it largely captures variance due to huge deletions

ncomp = 20
if(cancer == "LGG") ncomp <- 21
keepX = rep(10, ncomp)

log_tf_subset <- log(rnaseq_subset[match(TF_list, rownames(rnaseq_subset)),] + 1)

tf_sPC_list <- tf_loadings_list <- vector("list", length(TF_list) + 1)
names(tf_sPC_list) <- names(tf_loadings_list) <- c("all", TF_list)
for(tt in names(tf_sPC_list)) {
  cat("**************", tt)
  if(tt == "all") tf_dat <- log_tf_subset
  if(tt != "all") tf_dat <- log_tf_subset[-which(rownames(log_tf_subset) == tt),]
  tf_spca <- try(spca(t(tf_dat), ncomp = ncomp, center = TRUE, scale = TRUE,
    keepX = keepX))
  if(class(tf_spca)[1] == "try-error") {
    rr <- 16
    while(class(tf_spca)[1] == "try-error") {
      tf_spca <- try(spca(t(round(tf_dat, rr)), ncomp = ncomp, center = TRUE, 
        scale = TRUE, keepX = keepX))
      rr <- rr - 1
    }
    cat("***", tt, ": Rounded to", rr+1, "decimals\n")
  }
  tf_sPC <- tf_spca$x
  tf_loadings <- tf_spca$loading$X
  if(cancer == "LGG") {
    tf_sPC <- tf_sPC[,-1]
    tf_loadings <- tf_loadings[,-1]
  }
  keep_spca <- 1; index <- 2
  while(length(keep_spca) < number_sPCs) {
    tmp <- cor(tf_sPC[,c(keep_spca, index)])
    maxcor <- max(abs(tmp[upper.tri(tmp)]))
    ## Special case for TGCT and THYM, as we saw high correlation between
    ## TF variance components and methylation, so we don't include these TF sPCs
    if(cancer %in% c("XXXXXX")) {
      methyl_tmp <- apply(methyl_pergene_subset, 1, function(x) {
        methyl_tmp2 <- cor(cbind(x, tf_sPC[,c(index)]))[1,]
        return(methyl_tmp2[-1])
      })
      maxcor <- max(maxcor, max(abs(methyl_tmp), na.rm=TRUE))
    }
    if(maxcor > 0.3) index <- index + 1
    if(maxcor <= 0.3) {
      keep_spca <- c(keep_spca, index)
      index <- index + 1
      if(index > 20) break;
    }
  }
  tf_sPC_list[[tt]] <- tf_sPC[,keep_spca]
  tf_loadings_list[[tt]] <- tf_loadings[,keep_spca]
  cat(" :", keep_spca, "\n")
}
tf_keep_spca <- keep_spca + 1

#-----------------------------------------------------------------------------------
## Perform sPCA on log-transformed miRNA-seq data
## Choose first 5 uncorrelated sparse principal components
## (with other miRNA sPCs as well as the top 5 TF sPCs)

ncomp = 20
keepX = rep(10, ncomp)

log_mirna_subset <- log(mirna_subset + 1)
mirna_spca <- try(spca(t(log_mirna_subset), ncomp = ncomp, center = TRUE, scale = TRUE,
  keepX = keepX))
if(class(mirna_spca)[1] == "try-error") {
  rr <- 16
  while(class(mirna_spca)[1] == "try-error") {
    mirna_spca <- try(spca(t(round(log_mirna_subset, rr)), ncomp = ncomp, center = TRUE, 
      scale = TRUE, keepX = keepX))
    rr <- rr - 1
  }
  cat("*** mirna: Rounded to", rr+1, "decimals\n")
}
mirna_sPC <- mirna_spca$x
mirna_loadings <- mirna_spca$loading$X

## Remove any mir sPCs that are correlated with one another
keep_spca <- c(); index <- 1
other_sPCs <- tf_sPC_list
while(length(keep_spca) < number_sPCs) {
  tmp_other_sPCs <- lapply(other_sPCs, function(x) cbind(mirna_sPC[,index],x))
  tmp_cor <- lapply(tmp_other_sPCs, function(x) {
    xx <- cor(x)
    max(abs(xx[upper.tri(xx)]))
  })
  if(max(unlist(tmp_cor)) > 0.3) index <- index + 1
  if(max(unlist(tmp_cor)) <= 0.3) {
    keep_spca <- c(keep_spca, index)
    other_sPCs <- lapply(other_sPCs, function(x) cbind(x, mirna_sPC[,index]))
    index <- index + 1
    if(index > 20) stop("Not enough sPCs for miRNAs");
  }
}
mirna_sPC <- mirna_sPC[,keep_spca]
colnames(mirna_sPC) <- as.character(seq_len(number_sPCs))
mirna_loadings <- mirna_loadings[,keep_spca]
mir_keep_spca <- keep_spca


## Double-check correlation among all chosen sPCs
final_check <- lapply(tf_sPC_list, function(x) cbind(x[,1:5], mirna_sPC))
final_check_cor <- lapply(final_check, function(x) {
  xx <- cor(x)
  max(abs(xx[upper.tri(xx)]))
})
max(unlist(final_check_cor))



#-----------------------------------------------------------------------------------
## Fit regression to log(rna-seq + 1) values (already upper quartile normalized)
## for the first five principal components 

log_rnaseq_subset <- log(rnaseq_subset + 1)
rnaseq_pca <- pca(t(log_rnaseq_subset), ncomp = 5, center=TRUE, scale=TRUE)$x

rnaseq_resids0 <- matrix(NA, nrow=nrow(rnaseq_subset), ncol=ncol(rnaseq_subset))
rownames(rnaseq_resids0) <- rownames(rnaseq_subset)
colnames(rnaseq_resids0) <- substr(colnames(rnaseq_subset), 1, 12)
for(gg in seq_len(nrow(rnaseq_resids0))) {
  if(gg / 1000 == ceiling(gg / 1000)) cat(gg, "\n")
  lmdf <- data.frame(gene=log(rnaseq_subset[gg,]+1), rnaseq_pca)
  mod <- lm(gene ~ ., data=lmdf)
  rnaseq_resids0[gg,] <- residuals(mod)
}

#-----------------------------------------------------------------------------------
## Fit the multi-omic regression if needed
 
rnaseq_resids_fixedEffects <- matrix(NA, nrow=nrow(rnaseq_subset), ncol=ncol(rnaseq_subset))
rownames(rnaseq_resids_fixedEffects) <- rownames(rnaseq_subset)
colnames(rnaseq_resids_fixedEffects) <- substr(colnames(rnaseq_subset), 1, 12)

coef_fixedEffects <- matrix(NA, nrow=nrow(rnaseq_subset), ncol=(2*number_sPCs+3))
rownames(coef_fixedEffects) <- rownames(rnaseq_subset)
colnames(coef_fixedEffects) <- c("mut", "cna", "methyl", paste0("mirna.", seq_len(number_sPCs)),
  paste0("tf.", seq_len(number_sPCs)))
 
varComp_fixedEffects <- matrix(NA, nrow=nrow(rnaseq_subset), 
  ncol=(2*number_sPCs+4))
rownames(varComp_fixedEffects) <- rownames(rnaseq_subset)
colnames(varComp_fixedEffects) <- c("mut", "cna", "methyl", 
  paste0("mirna.", seq_len(number_sPCs)),
  paste0("tf.", seq_len(number_sPCs)), "residual")
 
for(gg in seq_len(nrow(rnaseq_resids0))) {
  if(gg / 1000 == ceiling(gg / 1000)) cat(gg, "\n")
  if(rownames(rnaseq_resids0)[gg] %in% TF_list) {
    tf_sPC <- tf_sPC_list[[rownames(rnaseq_resids0)[gg]]][,seq_len(number_sPCs)]
  }
  if(!rownames(rnaseq_resids0)[gg] %in% TF_list) {
    tf_sPC <- tf_sPC_list[["all"]][,seq_len(number_sPCs)]
  }
  colnames(tf_sPC) <- as.character(seq_len(number_sPCs))

  lmdf <- data.frame(rnaseq_resids0=rnaseq_resids0[gg,], mut=mut_pergene_subset[gg,],
    cna=cna_pergene_subset[gg,], methyl=methyl_pergene_subset[gg,],
    mirna=mirna_sPC[,seq_len(number_sPCs)], tf=tf_sPC)
  tmp <- paste(c(paste0("mirna.", seq_len(number_sPCs)),  
    paste0("tf.", seq_len(number_sPCs))), collapse=" + ")
  baseform <- paste0("rnaseq_resids0 ~ ", tmp)
  if(!sum(is.na(lmdf$mut))) baseform <- paste0(baseform, "+ mut")
  if(!sum(is.na(lmdf$cna))) baseform <- paste0(baseform, "+ cna")
  if(!sum(is.na(lmdf$methyl))) baseform <- paste0(baseform, "+ methyl")
  baseform <- as.formula(baseform) 
 
  mod <- lm(baseform, data=lmdf)

  ## Save coefficients and residuals
  cf <- coef(mod)
  index <- na.omit(match(names(cf), colnames(coef_fixedEffects)))
  coef_fixedEffects[gg, index] <- cf[-1]
  rnaseq_resids_fixedEffects[gg,] <- residuals(mod)

  ## Calculate variance components
  for(vc in names(cf)[-1]) {
    varComp_fixedEffects[gg, which(colnames(varComp_fixedEffects) == vc)] <- 
      var(cf[vc] * lmdf[,vc])
  }
  varComp_fixedEffects[gg, "residual"] <- var(mod$residuals)
}

#-----------------------------------------------------------------------------------
## Save results

if(subtype == "all") {  
  save(cancer, subtype, ga_keep_index, nmol, common_barcodes, mut_pergene_subset,
    clinical_subset, methyl_pergene_subset, cna_pergene_subset, rnaseq_subset, mirna_subset,
    TF_list, TF_full_list, file=paste0(basedir, cancer, "_results.RData"))
  save(rnaseq_resids0, rnaseq_resids_fixedEffects, coef_fixedEffects, varComp_fixedEffects,
    file=paste0(basedir, cancer, "_residuals.RData"))
  save(mirna_sPC, mirna_loadings, tf_sPC_list, tf_loadings_list, mir_keep_spca, tf_keep_spca,
    file=paste0(basedir, cancer, "_spca.RData"))
}


## BRCA subtype specific results
if(subtype != "all") {
  save(cancer, subtype, ga_keep_index, nmol, common_barcodes, mut_pergene_subset,
    clinical_subset, methyl_pergene_subset, cna_pergene_subset, rnaseq_subset, mirna_subset,
    TF_list, TF_full_list, file=paste0(basedir, cancer, "_", subtype, "_results.RData"))
  save(rnaseq_resids0, rnaseq_resids_fixedEffects, coef_fixedEffects, varComp_fixedEffects,
    file=paste0(basedir, cancer, "_", subtype, "_residuals.RData"))
  save(mirna_sPC, mirna_loadings, tf_sPC_list, tf_loadings_list, mir_keep_spca, tf_keep_spca,
    file=paste0(basedir, cancer, "_", subtype, "_spca.RData"))
}


