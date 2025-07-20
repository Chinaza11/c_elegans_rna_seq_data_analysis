#' ---
#' title: "C. elegans data analysis"
#' author: "Chinaza"
#' date: "`r Sys.Date()`"
#' ---

# =============================================================================
# ## setting working directory and loading packages
# =============================================================================
setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/differential_expression_analysis/deseq2/")

lst = c("dplyr", "ggplot2", "readr", "stringr", "RColorBrewer", "DESeq2", "EnhancedVolcano", "pheatmap")

load_libraries = function(list_of_packages){
  for (i in lst) {
    suppressPackageStartupMessages(library(i, character.only = TRUE))}
}

load_libraries(lst)

# =============================================================================
# ## setting my theme
# =============================================================================

my_theme = function() {
  base_size = 10
  line_colour = "grey85"
  text_colour = "grey20"
  theme_grey() %+replace%
    theme(
      line = element_line(colour = line_colour, linewidth = grid::unit(0.5, "pt"),
                          linetype = 1, lineend = "butt"),
      rect = element_rect(colour = NA, fill = NULL,
                          linewidth = grid::unit(0.5, "pt"), linetype = 1),
      text = element_text(colour = text_colour, face = "plain", size = 10,
                          lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0,
                          margin = margin()),
      axis.text = element_text(colour = text_colour, size = rel(0.8)),
      axis.ticks = element_line(colour = line_colour),
      legend.margin = margin(0, 0, 0, 0),
      legend.key = element_blank(),
      legend.key.size = unit(1.1 * base_size, "pt"),
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", colour = line_colour),
      panel.spacing.y = grid::unit(1, "lines"),
      plot.background = element_blank(),
      plot.tag = element_text(face = "bold", hjust = 0, vjust = 0),
      plot.tag.position = "topleft",
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0,
                                vjust = 0.5, margin = margin(b = base_size / 2)),
      strip.placement = "outside",
      strip.text = element_text(
        colour = text_colour,
        size = rel(0.8),
        margin = margin(base_size / 4, base_size / 4,
                        base_size / 4, base_size / 4, ))
    )
}

theme_set(my_theme())

# =============================================================================
# all functions
# =============================================================================

get_dds = function(metadata,counts){
  colData <- metadata
  countData <- counts %>% select(Geneid, all_of(colData$Sample))%>% filter(rowSums(pick(-Geneid)) >= 10) %>% data.frame()
  dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ Type, tidy = TRUE) %>% DESeq()
  return(dds)
}


get_pca = function(dds, title){
  pca <- dds %>%
    vst(blind = FALSE) %>%
    plotPCA(intgroup = "Type", returnData = TRUE) %>%
    mutate(name = str_remove(name, "^X")) %>%
    {ggplot(., aes(PC1, PC2, color = Type, label = name)) +
        geom_point() +
        coord_fixed() +
        scale_color_brewer(type = "qual", palette = "Dark2") +
        labs(
          title = title,
          x = { attr(., "percentVar")[1] * 100 } %>% round(digits = 1) %>% sprintf("PC1 (%.1f%%)", .),
          y = { attr(., "percentVar")[2] * 100 } %>% round(digits = 1) %>% sprintf("PC2 (%.1f%%)", .),
        )
    }
}


get_pca_with_samples_labelled = function(pca,title){
  pca + geom_label(show.legend = FALSE) + ggtitle(title)
}


getDeGenes_table = function(deseq_dataset, treatment1_name, treatment2_name){
  res = results(deseq_dataset, contrast=c("Type", treatment1_name, treatment2_name), alpha=0.05)
  print(paste("Result summary for", treatment1_name, "vs", treatment2_name))
  summary(res)
  #sum(res$padj < 0.05, na.rm=T)
  
  total = sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm=T)
  up = sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm=T)
  down = sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm=T)
  print(paste0("At padj < 0.05 and LFC > abs(1), total DEGs = ", total, " (",up," upregulated and ", down, " downregulated)")) 
  
  res_df = data.frame(res)
  return(res_df)
}


volcano_plot = function(dataframe, title){
  
  dataframe$padj = ifelse(-log10(dataframe$padj) < 20,dataframe$padj,1e-20)
  dataframe$log2FoldChange = ifelse(dataframe$log2FoldChange > 5, 5, dataframe$log2FoldChange)
  dataframe$log2FoldChange = ifelse(dataframe$log2FoldChange < -5, -5, dataframe$log2FoldChange)
  
  v = EnhancedVolcano(dataframe,
                      lab = rownames(dataframe),
                      x = 'log2FoldChange',
                      y = 'padj',
                      ylab = bquote(~-Log[10]~ italic(P-adj)),
                      pCutoff = 0.05,
                      FCcutoff = 1.0,
                      #pCutoffCol = "padj",
                      title = title,
                      subtitle = expression(italic('Horizontal dashed line cutoff at FDR 0.05 \nRight and left vertical dashed lines cutoff at LFC -1 and 1 respectively \nGenes with -Log10 P-adj value greater than 20 were set to 20 to allow better visualization')),
                      xlim = c(-6, 6),
                      ylim = c(0, 20),
                      legendLabels=c('NS', expression("Log"[2]*" FC"),'P-adj', expression('P-adj & Log'[2]*' FC')))
  return(v)
}


#function for descriptive label of the heatmap
add_suffix_to_label = function(sampleDistMatrix,pattern,suffix,row_or_col_names){
  if (row_or_col_names == "rownames"){
    for (i in 1:nrow(sampleDistMatrix)) {
      rowname <- rownames(sampleDistMatrix)[i]
      if (grepl(pattern, rowname)) {
        new_rowname <- paste0(rowname, suffix) 
        rownames(sampleDistMatrix)[i] <- new_rowname}}}
  
  if (row_or_col_names == "colnames"){
    for (i in 1:nrow(sampleDistMatrix)) {
      colname <- colnames(sampleDistMatrix)[i]
      if (grepl(pattern, colname)) {
        new_colname <- paste0(colname, suffix) 
        colnames(sampleDistMatrix)[i] <- new_colname}}}
  
  return(sampleDistMatrix)}


get_heatmap = function(dds,title){
  vsd = vst(dds, blind=FALSE)
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  rownames(sampleDistMatrix) = str_remove(rownames(sampleDistMatrix), "^X")
  colnames(sampleDistMatrix) = str_remove(colnames(sampleDistMatrix), "^X")
  
  patterns = list("1..1|1..2|1..3", "2..1|2..2|2..3", "3..1|3..2|3..3", "4..1|4..2|4..3")
  suffixes = list("_control", "_TYFD", "_TYEX", "_TPMG")
  
  for (i in 1:length(patterns)){
    sampleDistMatrix <- add_suffix_to_label(sampleDistMatrix,patterns[i],suffixes[i],"rownames")
    sampleDistMatrix <- add_suffix_to_label(sampleDistMatrix,patterns[i],suffixes[i],"colnames")}
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors, main=title)}


get_DEGs_barplot = function(TYFD,TYEX,TPMG,title){
  TYFD_up = sum(TYFD$padj < 0.05 & TYFD$log2FoldChange > 1, na.rm=T)
  TYFD_down = sum(TYFD$padj < 0.05 & TYFD$log2FoldChange < -1, na.rm=T)
  
  TYEX_up = sum(TYEX$padj < 0.05 & TYEX$log2FoldChange > 1, na.rm=T)
  TYEX_down = sum(TYEX$padj < 0.05 & TYEX$log2FoldChange < -1, na.rm=T)
  
  TPMG_up = sum(TPMG$padj < 0.05 & TPMG$log2FoldChange > 1, na.rm=T)
  TPMG_down = sum(TPMG$padj < 0.05 & TPMG$log2FoldChange < -1, na.rm=T)
  
  final_df = data.frame(regulation=rep(c("Upregulated","Downregulated"),each=3),
                        treatment=rep(c("TYFD","TYEX","TPMG"),2),
                        gene_count=c(TYFD_up,TYEX_up,TPMG_up,TYFD_down,TYEX_down,TPMG_down))
  
  ggplot(final_df, aes(x=treatment, y=gene_count, fill=regulation)) + 
    geom_bar(stat="identity", position=position_dodge2()) +
    geom_text(aes(label=gene_count), vjust = -0.2, position=position_dodge2(0.9)) +
    ggtitle(title)
}

# =============================================================================
# loading data
# =============================================================================

metadata <- read_tsv("metafile_edited.txt", col_types = "ccf") %>% relocate(Sample)

counts <- read_tsv("from_featurecounts/counts_edited.txt", col_types = "c_____") %>%
  dplyr::rename(all_of(metadata %>% select(Sample, File) %>% tibble::deframe()))

# =============================================================================
# DESeq2 analysis (all samples)
# =============================================================================

dds_all = get_dds(metadata,counts)

# =============================================================================
# PCA (all samples)
# =============================================================================

png("PCA.png", width=10, height=10, units="in", res=500)
pca_all = get_pca(dds_all,"C. elegans PCA, all samples")
pca_all
dev.off()

# =============================================================================
# Labelled PCA (all samples)
# =============================================================================

png("PCA_2.png", width=10, height=10, units="in", res=500)
get_pca_with_samples_labelled(pca_all,"C. elegans PCA, all samples with samples labelled")
dev.off()

# =============================================================================
# Heatmap (all samples)
# =============================================================================

png("heatmap.png", width=10, height=10, units="in", res=500)
get_heatmap(dds_all,"Heatmap, all samples")
dev.off()

# =============================================================================
# ## Control vs TPMG (all samples)
# =============================================================================
treatment1_name = "TPMG"
treatment2_name = "Control"

tpmg_vs_control_df = getDeGenes_table(dds_all, treatment1_name, treatment2_name)
#write.table(tpmg_vs_control_df, "tpmg_vs_control_df.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

tpmg_up_and_down = subset(tpmg_vs_control_df, padj < 0.05 & abs(log2FoldChange) > 1)
#write.table(tpmg_up_and_down, "tpmg_up_and_down.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("volcano_plot-TPMG.png", width=10, height=10, units="in", res=500)
tpmg_vs_control_plot = volcano_plot(tpmg_vs_control_df, paste(treatment1_name,"vs",treatment2_name, "\n"))
tpmg_vs_control_plot
dev.off()

# =============================================================================
# ## Control vs TYFD (all samples)
# =============================================================================
treatment1_name = "TYFD"
treatment2_name = "Control"

tyfd_vs_control_df = getDeGenes_table(dds_all, treatment1_name, treatment2_name)
#write.table(tyfd_vs_control_df, "tyfd_vs_control_df.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

tyfd_up_and_down = subset(tyfd_vs_control_df, padj < 0.05 & abs(log2FoldChange) > 1)
#write.table(tyfd_up_and_down, "tyfd_up_and_down.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("volcano_plot-TYFD.png", width=10, height=10, units="in", res=500)
tyfd_vs_control_plot = volcano_plot(tyfd_vs_control_df, paste(treatment1_name,"vs",treatment2_name, "\n"))
tyfd_vs_control_plot
dev.off()

# =============================================================================
# ## Control vs TYEX (all samples)
# =============================================================================
treatment1_name = "TYEX"
treatment2_name = "Control"

tyex_vs_control_df = getDeGenes_table(dds_all, treatment1_name, treatment2_name)
#write.table(tyex_vs_control_df, "tyex_vs_control_df.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

tyex_up_and_down = subset(tyex_vs_control_df, padj < 0.05 & abs(log2FoldChange) > 1)
#write.table(tyex_up_and_down, "tyex_up_and_down.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("volcano_plot-TYEX.png", width=10, height=10, units="in", res=500)
tyex_vs_control_plot = volcano_plot(tyex_vs_control_df, paste(treatment1_name,"vs",treatment2_name, "\n"))
tyex_vs_control_plot
dev.off()

# =============================================================================
# ## get barplot (all samples)
# =============================================================================

png("barplot.png", width=10, height=10, units="in", res=500)
get_DEGs_barplot(tyfd_vs_control_df,tyex_vs_control_df,tpmg_vs_control_df, "Differentially Expressed Genes (p-adj<0.05 & LFC>1), All Samples")
dev.off()


# =============================================================================
# =============================================================================
# ## sample 2..3 deleted
# =============================================================================
# =============================================================================

count_minus_2..3 = counts[, -which(names(counts) == "2--3")]
metadata_minus_2..3 = metadata[metadata$Sample != "2--3",]


dds_minus_2..3 = get_dds(metadata_minus_2..3,count_minus_2..3)

pca_minus_2..3 = get_pca(dds_minus_2..3,"C. elegans PCA, with sample 2..3 deleted")
pca_minus_2..3

get_pca_with_samples_labelled(pca_minus_2..3,"C. elegans PCA, samples labelled & sample 2..3 deleted")

get_heatmap(dds_minus_2..3,"Heatmap, with sample 2..3 deleted")


# ==========> TPMG <==========
treatment1_name = "TPMG"
treatment2_name = "Control"

tpmg_vs_control_df = getDeGenes_table(dds_minus_2..3, treatment1_name, treatment2_name)

tpmg_vs_control_plot = volcano_plot(tpmg_vs_control_df, paste(treatment1_name,"vs",treatment2_name, "\n"))
tpmg_vs_control_plot


# ==========> TYFD <==========
treatment1_name = "TYFD"
treatment2_name = "Control"

tyfd_vs_control_df = getDeGenes_table(dds_minus_2..3, treatment1_name, treatment2_name)

tyfd_vs_control_plot = volcano_plot(tyfd_vs_control_df, paste(treatment1_name,"vs",treatment2_name, "\n"))
tyfd_vs_control_plot


# ==========> TYEX <==========
treatment1_name = "TYEX"
treatment2_name = "Control"

tyex_vs_control_df = getDeGenes_table(dds_minus_2..3, treatment1_name, treatment2_name)

tyex_vs_control_plot = volcano_plot(tyex_vs_control_df, paste(treatment1_name,"vs",treatment2_name, "\n"))
tyex_vs_control_plot


# ==========> Barplot <==========
get_DEGs_barplot(tyfd_vs_control_df,tyex_vs_control_df,tpmg_vs_control_df, "Differentially Expressed Genes (p-adj<0.05 & LFC>1), with sample 2..3 deleted")

# =============================================================================
# 
# =============================================================================

sessionInfo()
