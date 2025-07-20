library(readr)
library(dplyr)
library(ggplot2)
library(DESeq2)

setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/differential_expression_analysis/deseq2/")

counts <- read_tsv("from_featurecounts/counts_edited.txt", col_types = "c_____")

metadata <- read_tsv("metafile_edited.txt", col_types = "ccf")

# =============================================================================
# function to get PCA plot for treatment vs control
# =============================================================================

doPCA <- function(counts, metadata, type) {
  colData <- metadata %>% filter(Type == "Control" | Type == type)
  countData  <- counts %>% select(Geneid, all_of(colData$File))%>% filter(rowSums(pick(-Geneid)) >= 10) %>% data.frame()
  
  output <- DESeqDataSetFromMatrix(countData, colData, design = ~ Type, tidy = TRUE)
  
  dds <- estimateSizeFactors(output)
  normcounts <- counts(dds, normalized = TRUE)
  
  deseq <- DESeq(output)
  
  vsd <- vst(deseq, blind = FALSE)
  pcaData <- plotPCA(vsd, intgroup = "Type", returnData = TRUE)
  
  plot <- ggplot(pcaData, aes(PC1, PC2, color = Type)) + 
    geom_point() + 
    coord_fixed() 
    
  return(plot)
}

theme_set(theme_classic() + theme(legend.title = element_blank()))

# =============================================================================
# TPMG
# =============================================================================

png("PCA-treatments_vs_control/TPMG.png", width=10, height=10, units="in", res=500)
plot1 <- doPCA(counts, metadata, "TPMG")
print(plot1)
dev.off()

# =============================================================================
# TYFD
# =============================================================================

png("PCA-treatments_vs_control/TYFD.png", width=10, height=10, units="in", res=500)
plot2 <- doPCA(counts, metadata, "TYFD")
print(plot2)
dev.off()

# =============================================================================
# TYEX
# =============================================================================

png("PCA-treatments_vs_control/TYEX.png", width=10, height=10, units="in", res=500)
plot3 <- doPCA(counts, metadata, "TYEX")
print(plot3)
dev.off()