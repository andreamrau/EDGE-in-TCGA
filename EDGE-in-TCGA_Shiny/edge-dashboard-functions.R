#---------------------------------------------------------------------------------------------
my_pheatmap <- function(gene_selection, cancer_selection, all_values, cancerabbrev, main="") {
  if(length(gene_selection) == 1) {
    sel <- all_values[gene %in% gene_selection & cancer %in% cancer_selection,
                      .(cancer, mut, cna, mirna, TF, methyl, genetic, residual)]
    sel <- as.data.frame(sel)
    rownames(sel) <- cancerabbrev$V2[match(sel$cancer, cancerabbrev$V1)]
    pve_values <- as.matrix(sel[,-1])

    annot <- all_values[gene %in% gene_selection & cancer %in% cancer_selection,
                        .(cancer, mean_rnaseq, mean_cna, mean_mut, mean_methyl, heritability)]
    annot <- as.data.frame(annot)
    rownames(annot) <- cancerabbrev$V2[match(annot$cancer, cancerabbrev$V1)]
  }
  if(length(gene_selection) > 1) {
    sel <- all_values[gene %in% gene_selection & cancer %in% cancer_selection,
                      .(gene, mut, cna, mirna, TF, methyl, genetic, residual)]
    sel <- as.data.frame(sel)
    rownames(sel) <- sel$gene
    pve_values <- as.matrix(sel[,-1])

    annot <- all_values[gene %in% gene_selection & cancer %in% cancer_selection,
                        .(gene, mean_rnaseq, mean_cna, mean_mut, mean_methyl, heritability)]
    annot <- as.data.frame(annot)
    rownames(annot) <- annot$gene
  }
  annot$mean_rnaseq <- log(annot$mean_rnaseq+1)
  annot <- annot[,-1]
  ann_colors = list(
    magma(21, begin=1, end=0),
    magma(21, begin=1, end=0),
    magma(21, begin=1, end=0),
    magma(21, begin=1, end=0),
    magma(21, begin=1, end=0)
  )
  
  colnames(annot) <- c("mean RNA-seq", "mean CNA", "mean mutation", "mean methylation", "heritability")
  names(ann_colors) <- colnames(annot)
  
  rmind <- c(which(colSums(is.na(annot))==nrow(annot)), which(apply(annot, 2, var) == 0))
  if(length(rmind)) annot <- annot[,-rmind]

  # ann_colors = list(
  #   c("white", "red"),
  #   c("white", "red"),
  #   c("#045A8D", "white", "red"),
  #   c("#045A8D", "white", "red"),
  #   c("white", "red")
  # )

  
  if(length(gene_selection) == 1) {
    pheatmap(pve_values, 
#           color = colorRampPalette(c("white", "red"))(25), 
             color = magma(21, begin=1, end=0), 
             cluster_cols=FALSE,
             display_numbers = matrix(ifelse(is.na(pve_values), "none", ""), nrow=nrow(pve_values),
                                     ncol=ncol(pve_values)),
             annotation_row=annot, 
            annotation_colors=ann_colors, 
            fontsize_row=8,
             labels_row = paste0(rownames(pve_values), "                 "),
             main=main)
    # ha <- rowAnnotation(annot,
    #                         col = list(mean_rnaseq = colorRamp2(c(0,8), c("white","red")),
    #                                    mean_cna = colorRamp2(c(-3,0,3), c("blue", "white","red")),
    #                                    mean_mut = colorRamp2(c(0,1), c("white", "red")),
    #                                    mean_methyl = colorRamp2(c(-3,0,3), c("blue", "white","red")),
    #                                    heritability = colorRamp2(c(0,1), c("white", "red"))),
    #                     width=unit(2, "cm"))
    # h <- Heatmap(pve_values, name="Value",
    #              col = colorRamp2(c(-1,seq(0,1,by=0.05)), c("grey", magma(21, begin=1, end=0))),
    #              show_row_names=TRUE,
    #              row_names_gp = gpar(fontsize = 6),
    #              heatmap_legend_param = list(at = c(-1, seq(0,1,by=0.1)),
    #                                          labels=c("none", as.character(seq(0,1,by=0.1)))))
    # draw(ha + h, row_dend_side="left")
    # decorate_annotation("mean_rnaseq", {grid.text("mean RNA-seq", unit(-10, "mm"), just="bottom", rot=90)})
                                                            
  }
  if(length(gene_selection) > 1) {
    pheatmap(pve_values, 
             color = magma(21, begin=1, end=0), 
  #           color = colorRampPalette(c("white", "red"))(25), 
             cluster_cols=FALSE,
             display_numbers = matrix(ifelse(is.na(pve_values), "none", ""), nrow=nrow(pve_values),
                                      ncol=ncol(pve_values)),
             annotation_row=annot, annotation_colors=ann_colors, fontsize_row=8,
             main=main)
  }

}

#---------------------------------------------------------------------------------------------
my_zoom <- function(gene_selection, cancer_selection, values, type="TF", horiz=FALSE) {
  if(values == "TF_values") {
    sqldb <- dbConnect(SQLite(), dbname = "dashboard_data/TF_values.sqlite")
    sel <- dbGetQuery(sqldb, paste0("select * from TF_values where gene ='", gene_selection, "' and cancer IN ",
                                    paste0("(", paste0("'", cancer_selection, "'", collapse=","), ")", collapse = "")))
    dbDisconnect(sqldb) ## Disconnect when done 
  }
  if(values == "mir_values") {
    sqldb <- dbConnect(SQLite(), dbname = "dashboard_data/mir_values.sqlite")
    sel <- dbGetQuery(sqldb, paste0("select * from mir_values where gene ='", gene_selection, "' and cancer IN ",
                                    paste0("(", paste0("'", cancer_selection, "'", collapse=","), ")", collapse = "")))
    dbDisconnect(sqldb) ## Disconnect when done 
  }
#  sel <- values[gene == gene_selection & cancer %in% cancer_selection]

  sel_list <- setNames(split(sel, seq(nrow(sel))), sel$cancer)
  sel_list <- lapply(sel_list, function(x) {
    if(type=="TF") {
      tmp <- x[1,3]
      colnames(x)[-c(1:3)] <- strsplit(TF_key[as.numeric(tmp),2], ";")$tfs
    }
    if(type=="mirs") {
      tmp <- x[1,2]
      colnames(x)[-c(1:3)] <- strsplit(mir_key[as.numeric(tmp),2], ";")$mirs
    }
    return(x)}
  )

  sel_list <- rbindlist(sel_list, idcol="cancer", fill=TRUE)
  sel_list <- data.frame(sel_list, stringsAsFactors=FALSE, check.names=FALSE)
  rownames(sel_list) <- sel_list$cancer
  sel_list <- sel_list[,-c(1:4)]
  sel_list[is.na(sel_list)] <- 0

  sel_list <- sel_list[,which(colSums(abs(sign(sel_list))) > 0)]
  ord <- apply(as.matrix(sel_list), 1, function(x) {
    tmp <- which(x!=0)
    if(type == "TF") return(data.frame(TF=names(x)[tmp], rank=rank(x[tmp], ties.method="random")))
    if(type == "mirs") return(data.frame(miRNA=names(x)[tmp], rank=rank(x[tmp], ties.method="random")))
  })
  ord <- rbindlist(ord, idcol="cancer")
  sel_list$cancer <- rownames(sel_list)
  if(type == "TF")  sel_list_tidy <- gather(sel_list, key=TF, value=loading, -cancer)
  if(type == "mirs") sel_list_tidy <- gather(sel_list, key=miRNA, value=loading, -cancer)
  sel_list_tidy <- sel_list_tidy[which(sel_list_tidy$loading != 0),]
  if(type == "TF")  sel_list_tidy <- left_join(sel_list_tidy, as.data.frame(ord), by=c("cancer", "TF"))
  if(type == "mirs") sel_list_tidy <- left_join(sel_list_tidy, as.data.frame(ord), by=c("cancer", "miRNA"))
  sel_list_tidy$loading2 <- sel_list_tidy$loading

  if(type == "TF" & horiz == FALSE) {
    p <- ggplot(sel_list_tidy) +
      geom_segment(aes(x=rank, xend=rank, y=0, yend=loading, color=loading2)) +
      geom_point(aes(x=rank, y=loading, color=loading2, label=TF)) +
      scale_colour_gradient2(low = "blue", high="red", mid="grey90", name="") +
      facet_wrap(~cancer, scales="free_x") + theme_bw() +
      geom_hline(yintercept=0) +
      labs(x=paste("TF Rank for", gene_selection), y="Weighted sPCA loading")
  }
  if(type == "mirs" & horiz == FALSE) {
    p <- ggplot(sel_list_tidy) +
      geom_segment(aes(x=rank, xend=rank, y=0, yend=loading, color=loading2)) +
      geom_point(aes(x=rank, y=loading, color=loading2, label=miRNA)) +
      scale_colour_gradient2(low = "blue", mid="grey90", high="red", name="") +
      facet_wrap(~cancer, scales="free_x") + theme_bw() +
      geom_hline(yintercept=0) +
      labs(x=paste("miRNA Rank for", gene_selection), y="Weighted sPCA loading")
  }
  if(type == "TF" & horiz == TRUE) {
    p <- ggplot(sel_list_tidy) +
      geom_segment(aes(y=rank, yend=rank, x=0, xend=loading, color=loading2)) +
      geom_point(aes(y=rank, x=loading, color=loading2, label=TF)) +
      geom_text(aes(y=rank, x=loading + 0.3*sign(loading), label=TF), size=2.5) +
      scale_colour_gradient2(low = "blue", high="red", mid="grey90", name="") +
      facet_wrap(~cancer, scales="free_x") + theme_minimal() +
      geom_vline(xintercept=0) +
      labs(y=" ", x="Weighted sPCA loading\n") + 
 #     xlim(-3, 3) +
      theme(axis.ticks = element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank()) 
  }
  if(type == "mirs" & horiz == TRUE) {
    p <- ggplot(sel_list_tidy) +
      geom_segment(aes(y=rank, yend=rank, x=0, xend=loading, color=loading2)) +
      geom_point(aes(y=rank, x=loading, color=loading2, label=miRNA)) +
      geom_text(aes(y=rank, x=loading + 0.3*sign(loading), label=miRNA), size=2.5) +
      scale_colour_gradient2(low = "blue", high="red", mid="grey90", name="") +
      facet_wrap(~cancer, scales="free_x") + theme_minimal() +
      geom_vline(xintercept=0) +
      labs(y=" ", x="Weighted sPCA loading\n") + 
#      xlim(-3, 3) +
      theme(axis.ticks = element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank())
  }
  return(p)
}

#---------------------------------------------------------------------------------------------
my_complexheatmap <- function(all_values, ind = NULL, 
                              molsource_choice = c("mut", "cna", "methyl", "mirna", "TF", "genetic"),
                              top200 = FALSE, rn=FALSE, top200_data=NULL) {
  
  if(top200 == FALSE) {
    ## Start with mutation
    heatmap_temp <- all_values[, c("gene", "cancer", molsource_choice[1]), with=FALSE]
    heatmap_temp <- heatmap_temp[which(heatmap_temp$gene %in% ind),] %>%
      mutate(cancer = as.factor(cancer)) %>%
      spread("cancer", molsource_choice[1])
    heatmap_temp$source <- molsource_choice[1]
    heatmap_dat <- heatmap_temp
    ## Add in other sources
    if(length(molsource_choice) > 1) {
      for(i in molsource_choice[-1]) {
        heatmap_temp <- all_values[, c("gene", "cancer", i), with=FALSE]
        heatmap_temp <- heatmap_temp[which(heatmap_temp$gene %in% ind),] %>%
          mutate(cancer = as.factor(cancer)) %>%
          spread("cancer", i)
        heatmap_temp$source <- i
        heatmap_dat <- rbind(heatmap_dat, heatmap_temp)
      }
    }
  }
  if(top200 == TRUE) {
  
    # heatmap_temp1 <- all_values %>% 
    #   dplyr::select(-chr, -strand, -tss, -residual) %>%
    #   gather(key=source, value=value, -gene, -cancer) 
    # choose_genes <- heatmap_temp1 %>%
    #   group_by(source, gene) %>%
    #   summarize(max_value = max(value, na.rm=TRUE)) %>%
    #   top_n(200) %>%
    #   ungroup() %>%
    #   dplyr::select(-max_value)
    # heatmap_dat <- left_join(choose_genes, heatmap_temp1, by=c("gene", "source")) %>%
    #   spread(key = cancer, value=value) %>%
    #   dplyr::select(gene, 3:19, source)
    # saveRDS(heatmap_dat, "dashboard_data/top200.rds")
  
    heatmap_dat <- top200_data %>%
      filter(source %in% molsource_choice)
  }

  heatmap_dat_orig <- heatmap_dat
  ## Replace NA's with -1 for now
  heatmap_dat[is.na(heatmap_dat)] <- -0.1
  ## Now do heatmap
  h <- Heatmap(heatmap_dat[,-c(1,ncol(heatmap_dat))], split=heatmap_dat$source, 
               name="Variance\ncomponent",
               col = colorRamp2(c(-0.1,0,seq(0.05,1,by=0.05)), c("grey", "white", magma(20, begin=1, end=0))),
              gap = unit(3, "mm"), show_row_names=FALSE,
               heatmap_legend_param = list(at = c(-0.1, seq(0,1,by=0.1)),
                                             labels=c("none", as.character(seq(0,1,by=0.1)))))
  if(rn ==TRUE) {
      ha <- rowAnnotation(gene = row_anno_text(heatmap_dat$gene,
                                               gp = gpar(fontsize = 7), offset = unit(0, "npc"),
                                               just="left"), width = max_text_width(heatmap_dat$gene))
      h <- h + ha
  }
  return(h)
}


#---------------------------------------------------------------------------------------------
my_newpheatmap <- function(gene_selection, cancer_selection, all_values, cancerabbrev, main="", rn=TRUE) {
  if(length(gene_selection) == 1) {
    sel <- all_values[gene %in% gene_selection & cancer %in% cancer_selection,
                      .(cancer, mut, cna, mirna, TF, methyl, genetic, residual)]
    sel <- as.data.frame(sel)
    rownames(sel) <- cancerabbrev$V2[match(sel$cancer, cancerabbrev$V1)]
    pve_values <- as.matrix(sel[,-1])
  }
  if(length(gene_selection) > 1) {
    sel <- all_values[gene %in% gene_selection & cancer %in% cancer_selection,
                      .(gene, mut, cna, mirna, TF, methyl, genetic, residual)]
    sel <- as.data.frame(sel)
    rownames(sel) <- sel$gene
    pve_values <- as.matrix(sel[,-1])
  }
  pve_values[is.na(pve_values)] <- -0.1
  colnames(pve_values) <- c("mutations", "copy number\nalterations", "miRNAs", "transcription\nfactors", "promoter\nmethylation",
                            "genetic", "residual")
  
    ## Replace NA's with -1 for now
    pve_values[is.na(pve_values)] <- -0.1
    ## Now do heatmap
    h <- Heatmap(pve_values,
                 name="Variance\ncomponent",
                 col = colorRamp2(c(-0.1,0,seq(0.05,1,by=0.05)), c("grey", "white", magma(20, begin=1, end=0))),
                 #           col = c("grey", "grey", magma(20, begin=1, end=0)),
                 #             na_col="white",
                 gap = unit(3, "mm"), show_row_names=FALSE,
                 heatmap_legend_param = list(at = c(-0.1, seq(0,1,by=0.1)),
                                             labels=c("none", as.character(seq(0,1,by=0.1)))))
    if(rn) {
      ha <- rowAnnotation(gene = row_anno_text(rownames(pve_values),
                                               offset = unit(0, "npc"),
                                               just="left"), width = max_text_width(rownames(pve_values)))
      h <- h + ha 
    }
      draw(h, heatmap_legend_side = "left")
}




#---------------------------------------------------------------------------------------------
my_dotplot <- function(all_values, gene_choice) {
  tmp <- all_values[, c("gene", "cancer", "mut", "cna", "methyl", "mirna", "TF", "genetic")]
  tmp <- tmp[which(tmp$gene == gene_choice),]
  tmp2 <- gather(tmp, key=type, value=value, -cancer, -gene)
  tmp2$cancer2 <- tmp2$cancer
  
  g <- ggplot(tmp2, aes(y=value, x=type, color=cancer2, fill=cancer2, label=cancer)) + 
    geom_segment(aes(y=0, 
                     yend=max(value, na.rm=TRUE)+.02, 
                     x=type, 
                     xend=type), color="grey",
                 linetype="dashed", 
                 size=0.25, show.legend=FALSE) +   # Draw dashed lines
    geom_jitter(size=4, show.legend=FALSE, height=0, width=0.2, alpha=0.65) + 
    geom_text(nudge_x = 0.35, size=3) +
    xlab("Molecular source of variation") +
    ylab("Variance component") + 
    scale_color_viridis(discrete=TRUE) +
    scale_fill_viridis(discrete=TRUE) +
    theme_classic() +
    guides(fill=FALSE, color=FALSE)  +
    ggtitle(gene_choice)
  return(g)
}
