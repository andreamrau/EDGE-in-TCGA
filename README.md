# EDGE-in-TCGA
Source code to reproduce results from ["Exploring Drivers of Gene Expression in The Cancer Genome Atlas"](https://www.biorxiv.org/content/early/2017/12/02/227926) by Rau et al. (2017)

This repository contains the following elements:

- **Annotation directory**: contains lists of transcription factors from the Ingenuity Pathway Analysis and TRRUST databases. Due to its size, the annotation file for the 450k human methylation array is not included here but can be found at https://support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html.

- **1_download_TCGA.R**: R script to download all available TCGA tumor samples from the Broad Firehose for RNA-seq (Illumina HiSeq), miRNA-seq (Illumina HiSeq and GA), somatic mutations, methylation (450K arrays), and copy number alterations data using the TCGA2STAT R package. Saves an .RData object containing all data in cancer-specific subdirectories.

- **2_format_and_preprocess_TCGA.R**: R script to read in the data downloaded from previous script, and perform pre-processing steps (match up 'omics data for each patient, select only Caucasian tumor samples, select most variable probes annotated with 1500bp of each gene's TSS, perform sparse principal component analysis on miRNA-seq data and expression data for TF's identifitied in IPA and TRRUST databases, fit initial per-gene model on RNA-seq data to regress out the first 5 principal components). Residuals from the pre-processing regression are saved in an .RData object in cancer-specific subdirectories.

- **3_GCTA_analysis_TCGA.R**: R script to fit a linear mixed model to standardized residuals from previous step using GCTA.

- **4_export_results_TCGA.R**: R script to read in results of previous analyses, calculate variance components for fixed effects, and export data for use in the Shiny app.

- **EDGE-in_TCGA_Shiny**: Directory containing all source code for creating the EDGE in TCGA Shiny app, found at http://andreamrau.shinyapps.io/edge-dashboard.

All of the R scripts above were run on a Slurm scheduler using the Rscript command, where the first argument represents the cancer type (BLCA, BRCA, CESC, ESCA, HNSC, KIRC, KIRP, LGG, LIHC, LUAD, PAAD, PCPG, PRAD, SARC, SKCM, STAD, THCA) and the second represents the cancer subtype (for all analyses here, we used "all"). Note that some minor modifications to the scripts (e.g., changing paths, creating cancer specific subdirectories, etc) are required prior to running the analyses.

If you use these scripts and/or our Shiny app in your research, please cite our paper:

Rau, A., Flister, M. J., Rui, H. and Livermore Auer, P. (2017) Exploring drivers of gene expression in The Cancer Genome Atlas. [bioRxiv](https://www.biorxiv.org/content/early/2017/12/02/227926), doi: https://doi.org/10.1101/227926.

This code, and the associated Shiny app, are distributed under the GNU public license (GPL-3).