#---------------------------------------------------------------------------------------------------
# Read in arguments from shell, load packages
args <- commandArgs(trailingOnly=TRUE)
cancer <- as.character(args[1])
# cancer <- "BLCA"
cat("***", cancer, "***\n")

user <- system("echo $USER", intern=TRUE)
if(user == "pauer") libpath <- "/raid-04/SPH/raua/R/x86_64-redhat-linux-gnu-library/3.2"
if(user == "raua")  libpath <- "/home/raua/Data/R/x86_64-redhat-linux-gnu-library/3.2"

library(TCGA2STAT, lib=libpath)

#---------------------------------------------------------------------------------------------------
# Reformat one TCGA2STAT command to allow for Illumina GA miRNA-seq data to be downloaded

my_miRNASeq <- function (ddoc, dlinks, dataset, platform, type = "count") 
{
    if (!(type %in% c("count", "rpmmm"))) {
        message("Error: Invalid type.")
        gdat <- NULL
        return(gdat)
    }
    if (type == "count") {
        type = "read_count"
    }
    if (type == "rpmmm") {
        type = "reads_per_million_miRNA_mapped"
    }
    keyWord = paste("", "Level_3__miR_gene_expression__data.Level_3", 
        sep = "")
    keyWord = paste("//a[contains(@href, '", keyWord, "')]", 
        sep = "")
    plinks = XML::xpathSApply(ddoc, keyWord, XML::xmlGetAttr, 
        "href")
    if(platform == "hiseq") {
      plinks = plinks[grepl(paste("*.", dataset, "[.]Merge_mirnaseq__.*.hiseq_mirnaseq__.*.tar[.]gz$", 
        sep = ""), plinks)]
    }
    if(platform == "ga") {
      plinks = plinks[grepl(paste("*.", dataset, "[.]Merge_mirnaseq__.*.ga_mirnaseq__.*.tar[.]gz$", 
        sep = ""), plinks)]
    }
    if (length(plinks) == 0) {
        message("Error: No data available for download. Please ensure the data is available from TCGA. \n")
        dat <- NULL
        return(dat)
    }
    timestamp <- unlist(strsplit(dlinks, "/"))
    timestamp <- timestamp[length(timestamp)]
    gdats <- list()
    for (i in 1:length(plinks)) {
        download_link = paste(dlinks, trim(plinks[i]), sep = "/")
        message("miRNAseq data will be imported! This may take some time!")
        utils::download.file(url = download_link, destfile = paste(dataset, 
            "-miRNAseqGene.tar.gz", sep = ""), method = "auto", 
            quiet = TRUE, mode = "w")
        fileList <- utils::untar(paste(dataset, "-miRNAseqGene.tar.gz", 
            sep = ""), list = TRUE)
        grepSearch = paste("*.", dataset, "[.]mirnaseq__.*.__Level_3__miR_gene_expression__data.data.txt$", 
            sep = "")
        fileList = fileList[grepl(grepSearch, fileList)]
        utils::untar(paste(dataset, "-miRNAseqGene.tar.gz", sep = ""), 
            files = fileList)
        fname = paste(dataset, "_", timestamp, "-miRNAseqGene.txt", 
            sep = "")
        file.rename(from = fileList, to = fname)
        file.remove(paste(dataset, "-miRNAseqGene.tar.gz", sep = ""))
        delFodler <- paste(getwd(), "/", strsplit(fileList, "/")[[1]][1], 
            sep = "")
        unlink(delFodler, recursive = TRUE)
        tmpCols = utils::read.delim(fname, nrows = 1, colClasses = "character")
        tmpdat = utils::read.delim(fname, skip = 1, sep = "\t", 
            stringsAsFactors = F)
        colOrder <- 1:ncol(tmpCols)
        colOrder <- colOrder[tmpCols[1, ] == type]
        gnames <- tmpdat[, 1]
        badg <- which(duplicated(gnames))
        if (length(badg) > 0) {
            gdat <- tmpdat[-badg, colOrder]
            colnames(gdat) <- colnames(tmpCols)[colOrder]
            colnames(gdat) <- gsub("\\.", "-", colnames(gdat))
            rownames(gdat) <- gnames[-badg]
        }
        if (length(badg) == 0) {
            gdat <- tmpdat[, colOrder]
            colnames(gdat) <- colnames(tmpCols)[colOrder]
            colnames(gdat) <- gsub("\\.", "-", colnames(gdat))
            rownames(gdat) <- gnames
        }
        message(paste(nrow(gdat), "genes have been imported!"))
        gdats[[i]] <- as.matrix(gdat)
        file.remove(fname)
    }
    if (length(gdats) == 1) {
        gdats <- gdats[[1]]
    }
    return(gdats)
}

my_getTCGA <- function (disease = "GBM", data.type = "RNASeq2", type = "", 
    filter = "Y", p = getOption("mc.cores", 2), clinical = FALSE, 
    cvars = "OS", platform="hiseq") 
{
    data.good <- c("RNASeq2", "RNASeq", "miRNASeq", "CNA_SNP", 
        "CNV_SNP", "CNA_CGH", "Methylation", "Mutation", "mRNA_Array", 
        "miRNA_Array")
    if (!(data.type %in% data.good)) {
        message("Error: Not recognized datatype for Firehose\n")
        dat <- NULL
        return(dat)
    }
    ldoc <- tryCatch({
        ldoc <- XML::htmlTreeParse("http://gdac.broadinstitute.org/runs/stddata__latest/", 
            useInternalNodes = T)
    }, error = function(e) {
        ldoc = NULL
    })
    if (is.null(ldoc)) {
        message("Error: Problem connect to Firehose. Please ensure Internet connection is working.\n")
        dat <- NULL
        return(dat)
    }
    datasets <- XML::xpathSApply(ldoc, "//a[contains(@href, 'Standardized+Data+Run+Release+Notes')]", 
        XML::xmlValue)
    if (!(disease %in% datasets)) {
        message("Error: Not recognized Disease Abbreviation for Firehose\n")
        dat <- NULL
        return(dat)
    }
    if (Sys.getenv("TAR") == "" | Sys.getenv("R_GZIPCMD") == 
        "") {
        message("Error: TAR is not installed in the system. Data unzip failed.\n")
        dat <- NULL
        return(dat)
    }
    dataset <- disease
    llinks = unlist(XML::xpathApply(ldoc, "//a[@href]", XML::xmlGetAttr, 
        "href"))
    dlinks = llinks[grepl(paste("/data/", dataset, "/", sep = ""), 
        llinks)]
    ddoc = XML::htmlTreeParse(dlinks, useInternalNodes = T)
    if (data.type == "miRNASeq") {
        if (type == "") {
            dat <- my_miRNASeq(ddoc = ddoc, dlinks = dlinks, dataset = dataset, platform = platform)
        }
        else {
            dat <- my_miRNASeq(ddoc = ddoc, dlinks = dlinks, dataset = dataset, 
                type = type, platform = platform)
        }
        if (is.null(dat)) {
            return(dat)
        }
        if (!clinical) {
            return(list(dat = dat, clinical = NULL, merged.dat = NULL))
        }
        if (clinical) {
            gdats <- list()
            cli <- Clinical(ddoc = ddoc, dlinks = dlinks, dataset = dataset)
            mdat <- NULL
            if (!is.null(cli)) {
                mdat <- MatrixMerge(dat = dat, cli = cli, cvars = cvars)
            }
            return(list(dat = dat, clinical = cli, merged.dat = mdat))
        }
    }
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)



#---------------------------------------------------------------------------------------------------
# Download data

## RSEM values for RNASeqV2
rnaseq <- getTCGA(disease=cancer, data.type="RNASeq2", clinical=TRUE)
clinical <- rnaseq$clinical
rnaseq <- rnaseq$dat

## By default, only Illumina HiSeq data are downloaded for V2 (RSEM values)
## We also add in Illumina GA data to increase sample size slightly
mirna_hiseq <- my_getTCGA(disease=cancer, data.type="miRNASeq", type="rpmmm", platform="hiseq")$dat
mirna_ga <- my_getTCGA(disease=cancer, data.type="miRNASeq", type="rpmmm", platform="ga")$dat

mut <- getTCGA(disease=cancer, data.type="Mutation", type="somatic")$dat
methyl <- getTCGA(disease=cancer, data.type="Methylation", type="450K")$dat
# Note that Y chromosomes are filtered out of CNA data
cna <- getTCGA(disease=cancer, data.type="CNA_SNP")$dat


#---------------------------------------------------------------------------------------------------
# Save data

remove <- ls()[which(!ls() %in% c("rnaseq", "clinical", "mirna_hiseq", "mirna_ga", 
  "mut", "methyl", "cna", "cancer"))]  
rm(list=remove)

save.image(paste0(cancer, "/TCGA_formattedData.RData"))

