args <- commandArgs(trailingOnly=TRUE)
cancer <- as.character(args[1])

### Load in Omics Data and Sparse PCAs
load(paste0("~/Data/TCGA_PROJECT/", cancer, "/genos_", cancer, ".RData"))
load(paste0("~/raua/TCGA_data/", cancer, "/", cancer, "_residuals.RData"))
load(paste0("~/raua/TCGA_data/", cancer, "/", cancer, "_spca.RData"))
load(paste0("~/raua/TCGA_data/", cancer, "/", cancer, "_results.RData"))
cancer <- as.character(args[1])

### Line up all the data 
sample.map <- subset(samples.keep, cel.codes %in% colnames(rnaseq_resids0) & filename %in% rownames(G))
G1 <- subset(G, rownames(G) %in% sample.map$filename)
sample.map <- sample.map[match(rownames(G1), sample.map$filename),]
rna_resids1 <- rnaseq_resids0[, match(sample.map$cel.codes, colnames(rnaseq_resids0))]
mut1 <- mut_pergene_subset[, match(sample.map$cel.codes, colnames(mut_pergene_subset))]

### Reformat CNA colnames
names.temp <- matrix(unlist(strsplit(colnames(cna_pergene_subset), split="-")), byrow=TRUE, ncol=7)
colnames(cna_pergene_subset) <- paste0(names.temp[,1], "-", names.temp[,2], "-", names.temp[,3])
cna1 <- cna_pergene_subset[, match(sample.map$cel.codes, colnames(cna_pergene_subset))]

### Reformat Methyl colnames
names.temp <- matrix(unlist(strsplit(colnames(methyl_pergene_subset), split="-")), byrow=TRUE, ncol=7)
colnames(methyl_pergene_subset) <- paste0(names.temp[,1], "-", names.temp[,2], "-", names.temp[,3])
meth1 <- methyl_pergene_subset[, match(sample.map$cel.codes, colnames(methyl_pergene_subset))]

### Reformat miRNA 
names.temp <- matrix(unlist(strsplit(rownames(mirna_sPC), split="-")), byrow=TRUE, ncol=8)
rownames(mirna_sPC) <- paste0(names.temp[,1], "-", names.temp[,2], "-", names.temp[,3])
mirna1 <- subset(mirna_sPC, rownames(mirna_sPC) %in% sample.map$cel.codes)
mirna2 <- mirna1[match(sample.map$cel.codes, rownames(mirna1)), ]

### Reformat TF_sPC
names.temp <- matrix(unlist(strsplit(rownames(tf_sPC_list[[1]]), split="-")), byrow=TRUE, ncol=7)
tf_list1 <- list()
for(i in 1:length(tf_sPC_list)){
rownames(tf_sPC_list[[i]]) <- paste0(names.temp[,1], "-", names.temp[,2], "-", names.temp[,3])
temp <- subset(tf_sPC_list[[i]], rownames(tf_sPC_list[[i]]) %in% sample.map$cel.codes)
tf_list1[[i]] <- temp[match(sample.map$cel.codes, rownames(temp)),]
}
names(tf_list1) <- names(tf_sPC_list)

### Read refGene file with all the TSS's for each gene
ref.gene <- read.table('~/Data/TCGA_PROJECT/refGene.txt', header=FALSE)

##################################################################
#### Andrea's final regression to obtain variance components #####
##################################################################
## Fit the multi-omic regression

coef_fixedEffects <- matrix(NA, nrow=nrow(rna_resids1), ncol=13)
rownames(coef_fixedEffects) <- rownames(rna_resids1)
colnames(coef_fixedEffects) <- c("mut", "cna", "methyl", paste0("mirna.", colnames(mirna2)[1:5]),
  paste0("tf.", 1:5))

varComp_Genetic <- matrix(NA, nrow=nrow(rna_resids1), ncol=3)
rownames(varComp_Genetic) <- rownames(rna_resids1)
colnames(varComp_Genetic) <- c("genetic", "residual", "total")

varComp_fixedEffects <- matrix(NA, nrow=nrow(rnaseq_subset), ncol=5)
rownames(varComp_fixedEffects) <- rownames(rnaseq_subset)
colnames(varComp_fixedEffects) <- c("mut", "cna", "methyl", "mirna", "tf")

for(gg in seq_len(nrow(rna_resids1))) {

  if(gg / 1000 == ceiling(gg / 1000)) cat(gg, "\n")

### Grab the correct TF sPCs
  if(rownames(rna_resids1)[gg] %in% TF_list) {
    tf_sPC <- tf_list1[[rownames(rna_resids1)[gg]]]
  }
  if(!rownames(rna_resids1)[gg] %in% TF_list) {
    tf_sPC <- tf_list1[["all"]]
  }
  colnames(tf_sPC) <- as.character(1:5)

### If there is only 1 carrier of a mutation, then change mut to all NAs 
  mut <- mut1[gg,]
  if(sum(mut==1, na.rm=TRUE)==1){
  mut <- rep(NA, times=length(mut))}

  resids.std <- rna_resids1[gg,]/sqrt(var(rna_resids1[gg,]))
### Construct the gene-level dataframe 
  lmdf <- data.frame(rna=resids.std, mut=mut, cna=cna1[gg,], methyl=meth1[gg,], mirna=mirna2[, 1:5], tf=tf_sPC[, 1:5])

### Run linear model without genetic effects
  tmp <- paste(paste(colnames(lmdf)[5:9], collapse=" + "), paste(colnames(lmdf)[10:14], collapse=" + "), sep=" + ")
  baseform <- paste0("rna ~ ", tmp)
  if(!sum(is.na(lmdf$mut))) baseform <- paste0(baseform, "+ mut")
  if(!sum(is.na(lmdf$cna))) baseform <- paste0(baseform, "+ cna")
  if(!sum(is.na(lmdf$methyl))) baseform <- paste0(baseform, "+ methyl")
  baseform <- as.formula(baseform)

  mod <- lm(baseform, data=lmdf)

################################################
### Run PLINK and GTCA with genetic effects ####
################################################

gene <- rownames(rna_resids1)[gg]
temp2 <- subset(ref.gene, ref.gene[,13] == gene)

### Only run this analysis if the gene is in the refgene file 
### Local GRM can only be run on autosomes
if(dim(temp2)[1]==0){

  ## Save coefficients
  cf <- coef(mod)
  index <- na.omit(match(names(cf), colnames(coef_fixedEffects)))
  coef_fixedEffects[gg, index] <- cf[-1]

  ## Calculate miRNA varComp
  vc <- names(cf)[grep(names(cf), pattern='mi')]
  varComp_fixedEffects[gg, "mirna"] <- var(as.matrix(lmdf[,vc]) %*% cf[vc])

  ### Calculate TF varComp
  vc <- names(cf)[grep(names(cf), pattern='tf')]
  varComp_fixedEffects[gg, "tf"] <- var(as.matrix(lmdf[,vc]) %*% cf[vc])

  ### Calculate any remaining varComps
  ## Calculate variance components
  names.remove <- c(names(cf)[grep(names(cf), pattern='tf')], names(cf)[grep(names(cf), pattern='mi')])
  names.keep <- setdiff(names(cf)[-1], names.remove)
  if(length(names.keep)>0){
  for(vc in names.keep) {
        varComp_fixedEffects[gg, which(colnames(varComp_fixedEffects) == vc)] <- var(cf[vc] * lmdf[,vc])
  }}

  varComp_Genetic[gg, "residual"] <- var(mod$residuals)
  varComp_Genetic[gg, "genetic"] <- 0
  varComp_Genetic[gg, "total"] <- var(mod$residuals)
}

if(dim(temp2)[1] > 0){

### Get all SNPs within 1MB on either side of the TSS.
cis.range <- c(temp2[,5] - 1000000, temp2[,5] + 1000000)
cis.snps <- subset(map, paste0("chr", map$chrom)==temp2[,3] & map$physical_pos >= cis.range[1] & map$physical_pos <= cis.range[2])

### If less than 30 markers or on chrX or chrY, don't estimate genetic effect
if(dim(cis.snps)[1] < 30 | !sum(temp2$V3 %in% paste0("chr", 1:22))){

  ## Save coefficients
  cf <- coef(mod)
  index <- na.omit(match(names(cf), colnames(coef_fixedEffects)))
  coef_fixedEffects[gg, index] <- cf[-1]

  ## Calculate miRNA varComp
  vc <- names(cf)[grep(names(cf), pattern='mi')]
  varComp_fixedEffects[gg, "mirna"] <- var(as.matrix(lmdf[,vc]) %*% cf[vc])

  ### Calculate TF varComp
  vc <- names(cf)[grep(names(cf), pattern='tf')]
  varComp_fixedEffects[gg, "tf"] <- var(as.matrix(lmdf[,vc]) %*% cf[vc])

  ### Calculate any remaining varComps
  ## Calculate variance components
  names.remove <- c(names(cf)[grep(names(cf), pattern='tf')], names(cf)[grep(names(cf), pattern='mi')])
  names.keep <- setdiff(names(cf)[-1], names.remove)
  if(length(names.keep)>0){
  for(vc in names.keep) {
        varComp_fixedEffects[gg, which(colnames(varComp_fixedEffects) == vc)] <- var(cf[vc] * lmdf[,vc])
  }}

  varComp_Genetic[gg, "residual"] <- var(mod$residuals)
  varComp_Genetic[gg, "genetic"] <- 0
  varComp_Genetic[gg, "total"] <- var(mod$residuals)
}else{

    ### Write out gene-specific map file
    write.map <- cbind(cis.snps[,c(2,4)], rep(0, times=dim(cis.snps)[1]), cis.snps[,3])
    write.table(write.map, paste0("~/Data/TCGA_PROJECT/", cancer, "/", gene, ".map"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

    ## Create .fam file and write out
    fam <- cbind(rownames(G1), rownames(G1), rep(0, times=dim(G1)[1]), rep(0, times=dim(G1)[1]), rep(2, times=dim(G1)[1]), rep(2, times=dim(G1)[1]))
    write.table(fam,  paste0("~/Data/TCGA_PROJECT/", cancer, "/", gene, ".fam"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

   ### Write out gene-specific .ped file
    dat2 <- G1[, which(colnames(G1) %in% cis.snps$man_fsetid)]
    dat3 <- cbind(fam, dat2)
    write.table(dat3, paste0("~/Data/TCGA_PROJECT/", cancer, "/", gene, ".ped"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")

    ### Run PLINK from the command line to convert to .bim and .bam for cis-analysis
    system(paste0("plink --file ~/Data/TCGA_PROJECT/", cancer, "/", gene, " --make-bed --noweb --out ~/Data/TCGA_PROJECT/", cancer, "/", gene))

    ### Build cis-GRM
    system(paste0("gcta64 --bfile ~/Data/TCGA_PROJECT/", cancer, "/", gene, " --make-grm --out ~/Data/TCGA_PROJECT/", cancer, "/", gene))

    ### Write out phenotype file
    phen.temp <- cbind(rownames(G1), rownames(G1), residuals(mod))
    write.table(phen.temp, paste0("~/Data/TCGA_PROJECT/", cancer, "/", gene, ".phen"), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

    ### Run cis-eQTL model 
    paste1 <- paste0("gcta64 --grm ~/Data/TCGA_PROJECT/", cancer, "/", gene, " --pheno ~/Data/TCGA_PROJECT/", cancer, "/", gene, ".phen --reml ")
    paste2 <- paste0("--out ~/Data/TCGA_PROJECT/", cancer, "/", gene, "_cis_out")
    system(paste(paste1, paste2))
    
    ### See if GCTA converged and wrote out results
    foo <- NULL
    foo <- system(paste0("ls ~/Data/TCGA_PROJECT/", cancer, "/ | grep  .hsq"), intern=TRUE)
    
    ### If REML didn't converge, then set genetic component to zero, and set residual and total variance to the variance of the residuals
    if(length(foo)==0){
    var.g <- matrix(0, nrow=3, ncol=3)
    var.g[,2] <- c(0, var(residuals(mod)), var(residuals(mod)))    
    }else{
    var.g <- read.table(paste0("~/Data/TCGA_PROJECT/", cancer, "/", gene, "_cis_out.hsq"), nrows=3, header=TRUE)
    }	  
    
    ### Delete all temporary gene specific files 
    system(paste0("rm ~/Data/TCGA_PROJECT/", cancer, "/", gene, "*"))

    ## Save coefficients and residuals
    cf <- coef(mod)
    index <- na.omit(match(names(cf), colnames(coef_fixedEffects)))
    coef_fixedEffects[gg, index] <- cf[-1]
    varComp_Genetic[gg,] <- var.g[,2]

    ### Calculate miRNA varComp
    vc <- names(cf)[grep(names(cf), pattern='mi')]
    varComp_fixedEffects[gg, "mirna"] <- var(as.matrix(lmdf[,vc]) %*% cf[vc])

    ### Calculate TF varComp	
    vc <- names(cf)[grep(names(cf), pattern='tf')]
    varComp_fixedEffects[gg, "tf"] <- var(as.matrix(lmdf[,vc]) %*% cf[vc])

    ### Calculate any remaining varComps
    ## Calculate variance components
    names.remove <- c(names(cf)[grep(names(cf), pattern='tf')], names(cf)[grep(names(cf), pattern='mi')])
    names.keep <- setdiff(names(cf)[-1], names.remove)
    if(length(names.keep)>0){
    for(vc in names.keep) {
    	   varComp_fixedEffects[gg, which(colnames(varComp_fixedEffects) == vc)] <- var(cf[vc] * lmdf[,vc])
    }}
}
}
cat("#########")
cat("\n")
cat(gg)
cat("\n")
cat(gene)
cat("\n")
cat("#########")
}

save(coef_fixedEffects, varComp_fixedEffects, varComp_Genetic, file=paste0("~/raua/TCGA_data/", cancer, "/", cancer, "_GCTA_results_standardized_July27.RData"))

