## Export data in appropriate format for Shiny viz

cancer_list <- c("BLCA", "BRCA", "CESC", "ESCA", "HNSC", "KIRC", 
  "KIRP", "LGG", "LIHC", "LUAD", "PAAD", "PCPG", "SARC", "SKCM", 
  "STAD", "THCA", "PRAD")

standardized <- TRUE
fullG <- FALSE

user <- system("echo $USER", intern=TRUE)
if(user == "pauer") {
  libpath <- "/raid-04/SPH/raua/R/x86_64-redhat-linux-gnu-library/3.2"
  basedir <- paste0("/raid-04/SPH/raua/TCGA_data/")
}
if(user == "raua")  {
  libpath <- "/home/raua/Data/R/x86_64-redhat-linux-gnu-library/3.2"
  basedir <- paste0("/home/raua/Data/TCGA_data/")
}
library(data.table, lib=libpath)
library(iotools, lib=libpath)

#--------------------------------------------------
## Clinical and sample size data
clinical <- vector("list", length(cancer_list))
names(clinical) <- cancer_list
nmol_all <- matrix(NA, nrow=7, ncol=length(cancer_list))
rownames(nmol_all) <- c("rnaseq", "mirna", "cna", "methyl", "mut", "samples", "mirna_ga")
colnames(nmol_all) <- cancer_list
for(cc in cancer_list) {
  cat("*** ", cc, "\n")
  load(paste0(basedir, cc, "/", cc, "_results.RData"))
  clinical[[cc]] <- clinical_subset
  nmol_all[-7,cc] <- as.numeric(nmol)
  if(!length(grep("_", cc))) nmol_all[7,cc] <- length(ga_keep_index)
}
saveRDS(clinical, file=paste0(basedir, "export_data/clinical.rds"), 
  compress=TRUE)
saveRDS(nmol_all, file=paste0(basedir, "export_data/nmol.rds"), 
  compress=TRUE)

#--------------------------------------------------
## Variance components
all <- vector("list", length(cancer_list))
names(all) <- cancer_list
for(cc in cancer_list) {
  cat("*** ", cc, "\n")
  load(paste0(basedir, cc, "/", cc, "_results.RData"))
  if(!standardized)
    load(paste0(basedir, cc, "/", cc, "_GCTA_results.RData"))
  if(standardized) {
    if(!fullG) load(paste0(basedir, cc, "/", cc, "_GCTA_results_standardized_July27.RData"))
    if(fullG) load(paste0(basedir, cc, "/", cc, "_GCTA_results_standardized_fullG_May25.RData"))
#   load(paste0(basedir, cc, "/", cc, "_residuals.RData"))
  }
  ## Check: since standardized residuals, the total variance should always < 1
  ## In some cases, there are really wacky values for the genetic and total components
  ## For these, we set the genetic component to 0
  if(standardized) {
    weird <- which(varComp_Genetic[,"total"] > 1)
    varComp_Genetic[weird,"genetic"] <- 0
#    varComp_Genetic[weird, "total"] <- varComp_Genetic[weird, "residual"]
  }
  all[[cc]] <- data.frame(gene=rownames(varComp_fixedEffects),
    mut=round(varComp_fixedEffects[,"mut"],3),
    cna=round(varComp_fixedEffects[,"cna"],3),
    methyl=round(varComp_fixedEffects[,"methyl"],3),
    mirna=round(varComp_fixedEffects[,"mirna"],3),
    TF=round(varComp_fixedEffects[,"tf"],3),
    genetic=round(varComp_Genetic[,"genetic"],3),
    residual=round(varComp_Genetic[,"residual"],3),
    heritability=round(varComp_Genetic[,"genetic"] / varComp_Genetic[,"total"],3),
    mean_rnaseq=round(rowMeans(rnaseq_subset),3),
    mean_cna=round(rowMeans(cna_pergene_subset),3),
    mean_mut=round(rowMeans(mut_pergene_subset),3),
    mean_methyl=round(rowMeans(methyl_pergene_subset),3), stringsAsFactors=FALSE)
  ## Now we rescale variance components to sum to 1
  all[[cc]][,2:8] <- all[[cc]][,2:8] / rowSums(all[[cc]][,2:8], na.rm=TRUE)
}
all_combine <- rbindlist(all, idcol="cancer")
if(!standardized) {
  saveRDS(all_combine, file=paste0(basedir, "export_data/varComp.rds"),
    compress=TRUE)
  write.csv.raw(all_combine, file=paste0(basedir, "export_data/varComp.csv")) 
}
if(standardized & !fullG) {
  saveRDS(all_combine, file=paste0(basedir, "export_data/varComp_standardized_July27.rds"),
    compress=TRUE)
  write.csv.raw(all_combine, file=paste0(basedir, "export_data/varComp_standardized_July27.csv")) 
}
if(standardized & fullG) {
  saveRDS(all_combine, file=paste0(basedir, "export_data/varComp_standardized_fullG_May25.rds"),
    compress=TRUE)
  write.csv.raw(all_combine, file=paste0(basedir, 
    "export_data/varComp_standardized_fullG_May25.csv")) 
}

#--------------------------------------------------
## miRNA sPCA results 

mirna_spca <- vector("list", length(cancer_list))
names(mirna_spca) <- cancer_list
mir_key <- matrix(NA, nrow=0, ncol=2)
colnames(mir_key) <- c("key", "mirs")
index <- 1
for(cc in cancer_list) {
  cat("*** ", cc, "\n")
  if(!standardized) 
    load(paste0(basedir, cc, "/", cc, "_GCTA_results.RData"))
  if(standardized) {
    if(!fullG) load(paste0(basedir, cc, "/", cc, "_GCTA_results_standardized_July27.RData"))
    if(fullG) load(paste0(basedir, cc, "/", cc, "_GCTA_results_standardized_fullG_May25.RData"))
  }
  load(paste0(basedir, cc, "/", cc, "_spca.RData"))

  mirna_coefs <- coef_fixedEffects[,grep("mirna", colnames(coef_fixedEffects))]
  mirna_loadings <- mirna_loadings[,1:5]
  mirna_loadings0 <- mirna_loadings[which(rowSums(sign(abs(mirna_loadings))) > 0),]

  ## We need a check in case the same mir appears in two components
  if(max(rowSums(abs(sign(mirna_loadings0)))) > 1) {
    double <- which(rowSums(abs(sign(mirna_loadings0)))>1)
    mirna_loadings_tmp <- mirna_loadings0[-double,]
    tmp <- mirna_loadings0[double,]
    if(is.vector(tmp)) {
      tmp <- matrix(tmp, nrow=1)
      rownames(tmp) <- names(double)
      colnames(tmp) <- colnames(mirna_loadings0)
    }
    tmp_bind <- do.call("rbind", lapply(data.frame(t(tmp)), function(x) {
      tmp1 <- which(x!=0) 
      tmp2 <- matrix(0, nrow=length(tmp1), ncol=length(x))
      tmp2[cbind(seq_len(length(tmp1)), tmp1)] <- x[tmp1]
      return(tmp2)
    }))
    colnames(tmp_bind) <- colnames(mirna_loadings_tmp)
    rownames(tmp_bind) <- rep(rownames(tmp), times=rowSums(abs(sign(tmp))))
    mirna_loadings0 <- rbind(mirna_loadings_tmp, tmp_bind)
  }
  if(nrow(mirna_loadings0) != 50) print("Whoa be careful ***********")

  tmp <- matrix(NA, nrow=nrow(coef_fixedEffects), ncol=nrow(mirna_loadings0))
  rownames(tmp) <- rownames(coef_fixedEffects) 
  colnames(tmp) <- rownames(mirna_loadings0)
  for(i in seq_len(nrow(tmp))) {
    tmp[i,] <- round(colSums(mirna_coefs[i,] * t(mirna_loadings0)), 3)
  }
  mirna_spca[[cc]] <- data.frame(gene=rownames(coef_fixedEffects), tmp,
    stringsAsFactors=FALSE, check.names=FALSE)
  mirs_in_pc <- paste0(colnames(mirna_spca[[cc]])[-1], collapse=";")
  key_exist <- grep(mirs_in_pc, mir_key[,"mirs"], fixed=TRUE)
  if(!length(key_exist)) {
    mirna_spca[[cc]] <- data.frame(key=rep(index, nrow(mirna_spca[[cc]])), 
      mirna_spca[[cc]])
    mir_key <- rbind(mir_key, c(index, mirs_in_pc))
    index <- index + 1
  }
  if(length(key_exist) == 1) {
    mirna_spca[[cc]] <- data.frame(key=rep(key_exist, nrow(mirna_spca[[cc]])),
      mirna_spca[[cc]])
  }
  if(length(key_exist) > 1) print("We have a problem")
  colnames(mirna_spca[[cc]])[-c(1:2)] <- paste0("mir", seq_len(ncol(mirna_spca[[cc]])-2))
}
mirna_spca_combine <- rbindlist(mirna_spca, idcol="cancer", fill=TRUE)
mirna_spca_combine[is.na(mirna_spca_combine)] <- 0
if(!standardized) {
  saveRDS(mirna_spca_combine, file=paste0(basedir, "export_data/mirna_spca.rds"),
    compress=TRUE)
  saveRDS(mir_key, file=paste0(basedir, "export_data/mirna_spca_key.rds"),
    compress=TRUE)
  write.csv.raw(mir_key, file=paste0(basedir, "export_data/mirna_spca_key.csv"))
}
if(standardized & !fullG) {
  saveRDS(mirna_spca_combine, file=paste0(basedir, 
    "export_data/mirna_spca_standardized_July27.rds"),
    compress=TRUE)
  saveRDS(mir_key, file=paste0(basedir, 
    "export_data/mirna_spca_key_standardized_July27.rds"),
    compress=TRUE)
  write.csv.raw(mir_key, file=paste0(basedir, 
    "export_data/mirna_spca_key_standardized_July27.csv"))
}
if(standardized & fullG) {
  saveRDS(mirna_spca_combine, file=paste0(basedir, 
    "export_data/mirna_spca_standardized_fullG_May25.rds"),
    compress=TRUE)
  saveRDS(mir_key, file=paste0(basedir, 
    "export_data/mirna_spca_key_standardized_fullG_May25.rds"),
    compress=TRUE)
  write.csv.raw(mir_key, file=paste0(basedir, 
    "export_data/mirna_spca_key_standardized_fullG_May25.csv"))
}

(apply(mir_key, 1, function(x) length(unlist(strsplit(x[2], split=";")))))
table(apply(mir_key, 1, function(x) length(unlist(strsplit(x[2], split=";")))))



#--------------------------------------------------
## TF sPCA results

tf_spca <- vector("list", length(cancer_list))
names(tf_spca) <- cancer_list
tf_key <- matrix(NA, nrow=0, ncol=2)
colnames(tf_key) <- c("key", "tfs")
index <- 1
for(cc in cancer_list) {
  cat("*** ", cc, "\n")
  if(!standardized)
    load(paste0(basedir, cc, "/", cc, "_GCTA_results.RData"))
  if(standardized) {
    if(!fullG) load(paste0(basedir, cc, "/", cc, "_GCTA_results_standardized_July27.RData"))
    if(fullG) load(paste0(basedir, cc, "/", cc, "_GCTA_results_standardized_fullG_May25.RData"))
  }
  load(paste0(basedir, cc, "/", cc, "_spca.RData"))
  tf_coefs <- coef_fixedEffects[,grep("tf.", colnames(coef_fixedEffects))]
  tf_list <- names(tf_loadings_list)[-1]

  tmp <- vector("list", nrow(tf_coefs))
  names(tmp) <- rownames(tf_coefs)
  for(g in seq_len(nrow(tf_coefs))) {
    gene <- rownames(tf_coefs)[g]
    if(gene %in% tf_list) {
      tf_loadings <- tf_loadings_list[[gene]][,1:5]
    }
    if(!gene %in% tf_list) {
      tf_loadings <- tf_loadings_list[["all"]][,1:5]
    }
    tf_loadings0 <- tf_loadings[which(rowSums(sign(abs(tf_loadings))) > 0),]

    ## We need a check in case the same TF appears in two components
    if(max(rowSums(abs(sign(tf_loadings0)))) > 1) {
      double <- which(rowSums(abs(sign(tf_loadings0)))>1)
      tf_loadings_tmp <- tf_loadings0[-double,]
      tmpp <- tf_loadings0[double,]
      if(is.vector(tmpp)) {
        tmpp <- matrix(tmpp, nrow=1)
        rownames(tmpp) <- names(double)
        colnames(tmpp) <- colnames(tf_loadings0)
      }
      tmp_bind <- do.call("rbind", lapply(data.frame(t(tmpp)), function(x) {
        tmp1 <- which(x!=0) 
        tmp2 <- matrix(0, nrow=length(tmp1), ncol=length(x))
        tmp2[cbind(seq_len(length(tmp1)), tmp1)] <- x[tmp1]
        return(tmp2)
      }))
      colnames(tmp_bind) <- colnames(tf_loadings_tmp)
      rownames(tmp_bind) <- rep(rownames(tmpp), times=rowSums(abs(sign(tmpp))))
      tf_loadings0 <- rbind(tf_loadings_tmp, tmp_bind)
    }
    if(nrow(tf_loadings0) != 50) print("Whoa be careful ***********")

    tmp[[g]] <- matrix(NA, nrow=1, ncol=nrow(tf_loadings0))
    rownames(tmp[[g]]) <- gene 
    colnames(tmp[[g]]) <- rownames(tf_loadings0)
    tmp[[g]][1,] <- round(colSums(tf_coefs[g,] * t(tf_loadings0)), 3) 
    tmp[[g]] <- as.data.frame(tmp[[g]], check.names=FALSE, stringsAsFactors=FALSE)

    tfs_in_pc <- paste0(colnames(tmp[[g]]), collapse=";")
    key_exist <- grep(tfs_in_pc, tf_key[,"tfs"], fixed=TRUE)
    if(!length(key_exist)) {
      tmp[[g]] <- data.frame(key=rep(index, nrow(tmp[[g]])), tmp[[g]])
      tf_key <- rbind(tf_key, c(index, tfs_in_pc))
      index <- index + 1
    }
    if(length(key_exist) == 1) {
      tmp[[g]] <- data.frame(key=rep(key_exist, nrow(tmp[[g]])), tmp[[g]])
    }
    if(length(key_exist) > 1) print("We have a problem")
    colnames(tmp[[g]])[-1] <- paste0("tf", seq_len(ncol(tmp[[g]])-1))
  }
  tf_spca[[cc]] <- rbindlist(tmp, fill=TRUE, idcol="gene")
}
tf_spca_combine <- rbindlist(tf_spca, idcol="cancer", fill=TRUE)
tf_spca_combine[is.na(tf_spca_combine)] <- 0
if(!standardized) {
  saveRDS(tf_spca_combine, file=paste0(basedir, "export_data/tf_spca.rds"),
    compress=TRUE)
  write.csv.raw(tf_spca_combine, file=paste0(basedir, "export_data/tf_spca.csv"))
  saveRDS(tf_key, file=paste0(basedir, "export_data/tf_spca_key.rds"),
    compress=TRUE)
  write.csv.raw(tf_key, file=paste0(basedir, "export_data/tf_spca_key.csv"))
}
if(standardized & !fullG) {
  saveRDS(tf_spca_combine, file=paste0(basedir, "export_data/tf_spca_standardized_July27.rds"),
    compress=TRUE)
  write.csv.raw(tf_spca_combine, file=paste0(basedir, "export_data/tf_spca_standardized_July27.csv"))
  saveRDS(tf_key, file=paste0(basedir, "export_data/tf_spca_key_standardized_July27.rds"),
    compress=TRUE)
  write.csv.raw(tf_key, file=paste0(basedir, "export_data/tf_spca_key_standardized_July27.csv"))
}
if(standardized & fullG) {
  saveRDS(tf_spca_combine, file=paste0(basedir, 
    "export_data/tf_spca_standardized_fullG_May25.rds"),
    compress=TRUE)
  write.csv.raw(tf_spca_combine, file=paste0(basedir, 
    "export_data/tf_spca_standardized_fullG_May25.csv"))
  saveRDS(tf_key, file=paste0(basedir, 
    "export_data/tf_spca_key_standardized_fullG_May25.rds"),
    compress=TRUE)
  write.csv.raw(tf_key, file=paste0(basedir, 
    "export_data/tf_spca_key_standardized_fullG_May25.csv"))
}


table(apply(tf_key, 1, function(x) length(unlist(strsplit(x[2], split=";")))))


#--------------------------------------------------
## cor list of components across cancers
library(data.table)
library(tidyr)
library(dplyr)

all_values <- readRDS(paste0(basedir, "export_data/varComp_standardized_July27.rds"))
cancer_choice <- unique(all_values$cancer)
gene_names <- unique(all_values$gene)
varcomp <- all_values %>% select(-starts_with("mean"), -heritability)
varcomp_long <- gather(varcomp, source, varcomp, -cancer, -gene)
varcomp_long[is.na(varcomp_long)] <- 0
varcomp_long$source <- factor(varcomp_long$source, levels=colnames(varcomp)[-c(1:2)])

varcomp_DT <- setDT(varcomp_long)
cor_list <- vector("list", 7)
names(cor_list) <- colnames(varcomp)[-c(1:2)]
for(vc in colnames(varcomp)[-c(1:2)]) {
  cat(vc, "\n")
  tmp <- varcomp_DT[source == vc, .(cancer, gene, varcomp)]
  tmp_wide <- dcast.data.table(tmp, gene~cancer) %>%
    na.omit()
  tmp_wide1 <- tmp_wide
  tmp_wide1$gene <- NULL

  if(vc == "genetic") {
    index <- which(colSums(tmp_wide1) == 0)
    if(length(index)) tmp_wide1 <- tmp_wide1[,-index]
  }
  corrcalc <- cor(tmp_wide1)
  cor_list[[vc]] <- corrcalc
}
names(cor_list) <- c("mutation", "CNA", "methylation", "miRNA", "TF", "genetic", "residual")
saveRDS(file=paste0(basedir, "export_data/corlist_July27.rds"), cor_list)

#--------------------------------------------------
## cor list of pairwise variance components for each cancer
cor_list2 <- matrix(NA, nrow = 15, ncol = length(unique(varcomp$cancer)))
colnames(cor_list2) <- unique(varcomp$cancer)
source_choice <- c("mut", "cna", "methyl", "mirna", "TF", "genetic")
nams <- c()
for(s1 in 1:5) {
  for(s2 in (s1+1):6) {
    source1 <- source_choice[s1]
    source2 <- source_choice[s2]
    nams <- c(nams, paste0(source1, ".", source2))
  }
}
rownames(cor_list2) <- nams
for(can in unique(varcomp$cancer)) {
  cat(can, "\n")
  tmp <- varcomp_DT[cancer == can, .(gene, source, varcomp)]
  tmp_wide <- dcast.data.table(tmp, gene~source) %>%
    na.omit()
  for(s1 in 1:5) {
    for(s2 in (s1+1):6) {
      source1 <- source_choice[s1]
      source2 <- source_choice[s2]
      tmp_wide_choice <- tmp_wide %>%
        select_(source1, source2) %>%
        as.matrix
      corrcalc <- cor(tmp_wide_choice)[1,2]
      rowindex <- which(rownames(cor_list2) == paste0(source1, ".", source2))
      cor_list2[rowindex, which(colnames(cor_list2) == can)] <- corrcalc
    }
  }
}
saveRDS(file=paste0(basedir, "export_data/corlist2_July27.rds"), cor_list2)















## Save a copy of old results (very sparse matrices)

#--------------------------------------------------
## miRNA sPCA results 

mirna_spca <- vector("list", length(cancer_list))
names(mirna_spca) <- cancer_list
for(cc in cancer_list) {
  cat("*** ", cc, "\n")
  load(paste0(basedir, cc, "/", cc, "_residuals.RData"))
  load(paste0(basedir, cc, "/", cc, "_spca.RData"))
  mirna_coefs <- coef_fixedEffects[,4:8]

  mirna_loadings0 <- mirna_loadings[which(rowSums(sign(abs(mirna_loadings))) > 0),]
  tmp <- matrix(NA, nrow=nrow(coef_fixedEffects), ncol=nrow(mirna_loadings0))
  rownames(tmp) <- rownames(coef_fixedEffects) 
  colnames(tmp) <- rownames(mirna_loadings0)
  for(i in seq_len(nrow(tmp))) {
    tmp[i,] <- round(colSums(mirna_coefs[i,] * t(mirna_loadings0)), 1)
  }
  mirna_spca[[cc]] <- data.frame(gene=rownames(coef_fixedEffects), tmp,
    stringsAsFactors=FALSE, check.names=FALSE)
}
mirna_spca_combine <- rbindlist(mirna_spca, idcol="cancer", fill=TRUE)
mirna_spca_combine[is.na(mirna_spca_combine)] <- 0
saveRDS(mirna_spca_combine, file=paste0(basedir, "export_data/mirna_spca.rds"),
  compress=TRUE)
write.csv.raw(mirna_spca_combine, file=paste0(basedir, "export_data/mirna_spca.csv"))

## TODO: Subset to remove elements with max|x| < some threshold

#--------------------------------------------------
## TF sPCA results

tf_spca <- vector("list", length(cancer_list))
names(tf_spca) <- cancer_list
for(cc in cancer_list) {
  cat("*** ", cc, "\n")
  load(paste0(basedir, cc, "/", cc, "_residuals.RData"))
  load(paste0(basedir, cc, "/", cc, "_spca.RData"))
  tf_coefs <- coef_fixedEffects[,9:13]
  tf_list <- names(tf_loadings_list)[-1]

  tmp <- vector("list", nrow(tf_coefs))
  names(tmp) <- rownames(tf_coefs)
  for(g in seq_len(nrow(tf_coefs))) {
    gene <- rownames(tf_coefs)[g]
    if(gene %in% tf_list) {
      tf_loadings <- tf_loadings_list[[gene]]
    }
    if(!gene %in% tf_list) {
      tf_loadings <- tf_loadings_list[["all"]]
    }
    tf_loadings0 <- tf_loadings[which(rowSums(sign(abs(tf_loadings))) > 0),]
    tmp[[g]] <- matrix(NA, nrow=1, ncol=nrow(tf_loadings0))
    rownames(tmp[[g]]) <- gene 
    colnames(tmp[[g]]) <- rownames(tf_loadings0)
    tmp[[g]][1,] <- round(colSums(tf_coefs[g,] * t(tf_loadings0)), 1) 
    tmp[[g]] <- as.data.frame(tmp[[g]], check.names=FALSE, stringsAsFactors=FALSE)
  }
  tf_spca[[cc]] <- rbindlist(tmp, fill=TRUE, idcol="gene")
}
tf_spca_combine <- rbindlist(tf_spca, idcol="cancer", fill=TRUE)
tf_spca_combine[is.na(tf_spca_combine)] <- 0
saveRDS(tf_spca_combine, file=paste0(basedir, "export_data/tf_spca.rds"),
  compress=TRUE)
#write.csv.raw(tf_spca_combine, file=paste0(basedir, "export_data/tf_spca.csv"))

## TODO: Subset to remove elements with max|x| < some threshold