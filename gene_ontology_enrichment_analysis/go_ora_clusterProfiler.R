# =============================================================================
# ## setting working directory and loading packages
# =============================================================================
setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result")

library(clusterProfiler)
library(org.Ce.eg.db)
library(tibble)
library(enrichplot)
library(ggplot2)


go_analysis = function(path) {

  df = read.table(path, header=TRUE, sep="\t")
  df = rownames_to_column(df, var = "ncbi_locus_tag")
  #print(head(df))
  
  ids = read.csv("DAVID-locus_tag_to_entrez_id.csv")
  names(ids)[names(ids) == "From"] = "ncbi_locus_tag"
  names(ids)[names(ids) == "To"] = "Entrez_id"
  #print(head(ids))
  
  merged_df = merge(df, ids, by = "ncbi_locus_tag", all = TRUE)
  #print(head(merged_df))
  
  genes_of_interest = as.character((subset(merged_df, padj < 0.05 & abs(log2FoldChange) > 1))$Entrez_id)
  
  ego = enrichGO(gene          = genes_of_interest,
                  universe      = as.character(merged_df$Entrez_id),
                  OrgDb         = org.Ce.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  result = data.frame(ego)
  #print(head(result))
  
  modify_plot = theme(plot.title = element_text(size = 24, face ="bold", hjust = 0.5))
  
  dotplot1 = dotplot(ego, showCategory = dim(result)[1]) + 
    labs(title="Enriched GO Terms for Biological Process (all terms)", 
         caption="GeneRatio ==> ratio of differential expressed genes that are annotated in a term\nCount ==> GeneCount") + 
    modify_plot
  
  ## use simplify to remove redundant terms
  ego2 = simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  result2 = data.frame(ego2)
  #print(head(result2))
  
  dotplot2 = dotplot(ego2, showCategory = dim(result2)[1]) + 
          labs(title="Enriched GO Terms for Biological Process (redundant terms removed)",
               caption="GeneRatio ==> ratio of differential expressed genes that are annotated in a term\nCount ==> GeneCount") + 
          modify_plot 
  
  go_plot = goplot(ego2, showCategory = dim(result2)[1]) + 
    ggtitle("Directed Acyclic Graph of Enrichment Analysis") + 
    modify_plot
  
  cnet_plot = cnetplot(ego2, showCategory = dim(result2)[1]) + 
    ggtitle("Network Plot of Enriched Terms") + 
    modify_plot
  
  ego2 = pairwise_termsim(ego2) 
  emap_plot = emapplot(ego2) + 
    ggtitle("Enrichment Map showing overlapping terms") + 
    modify_plot
  
  output = list(result_all_terms = result, dotplot_all_terms = dotplot1, result_redundant_terms_removed = result2, dotplot_redundant_terms_removed = dotplot2, go_plot = go_plot, cnet_plot = cnet_plot, emap_plot = emap_plot)
  
  return(output)
}

# =============================================================================
# ## Control vs TPMG 
# =============================================================================

path = "differential_expression_analysis/deseq2/tpmg_vs_control_df.txt"

TPMG = go_analysis(path)


#write.table(TPMG$result_all_terms, "gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TPMG/result_all_terms.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(TPMG$result_redundant_terms_removed, "gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TPMG/result_redundant_terms_removed.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TPMG/dotplot_all_terms.png", width=14, height=14, units="in", res=500)
print(TPMG$dotplot_all_terms)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TPMG/dotplot_redundant_terms_removed.png", width=14, height=14, units="in", res=500)
print(TPMG$dotplot_redundant_terms_removed)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TPMG/go_plot.png", width=10, height=10, units="in", res=500)
print(TPMG$go_plot)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TPMG/cnet_plot.png", width=10, height=10, units="in", res=500)
print(TPMG$cnet_plot)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TPMG/emap_plot.png", width=10, height=10, units="in", res=500)
print(TPMG$emap_plot)
dev.off()

# =============================================================================
# ## Control vs TYFD
# =============================================================================

path = "differential_expression_analysis/deseq2/tyfd_vs_control_df.txt"

TYFD = go_analysis(path)

#write.table(TYFD$result_all_terms, "gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYFD/result_all_terms.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(TYFD$result_redundant_terms_removed, "gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYFD/result_redundant_terms_removed.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYFD/dotplot_all_terms.png", width=14, height=14, units="in", res=500)
print(TYFD$dotplot_all_terms)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYFD/dotplot_redundant_terms_removed.png", width=14, height=14, units="in", res=500)
print(TYFD$dotplot_redundant_terms_removed)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYFD/go_plot.png", width=10, height=10, units="in", res=500)
print(TYFD$go_plot)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYFD/cnet_plot.png", width=10, height=10, units="in", res=500)
print(TYFD$cnet_plot)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYFD/emap_plot.png", width=10, height=10, units="in", res=500)
print(TYFD$emap_plot)
dev.off()

# =============================================================================
# ## Control vs TYEX
# =============================================================================

path = "differential_expression_analysis/deseq2/tyex_vs_control_df.txt"

TYEX = go_analysis(path)

#write.table(TYEX$result_all_terms, "gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYEX/result_all_terms.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(TYEX$result_redundant_terms_removed, "gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYEX/result_redundant_terms_removed.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYEX/dotplot_all_terms.png", width=14, height=14, units="in", res=500)
print(TYEX$dotplot_all_terms)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYEX/dotplot_redundant_terms_removed.png", width=14, height=14, units="in", res=500)
print(TYEX$dotplot_redundant_terms_removed)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYEX/go_plot.png", width=10, height=10, units="in", res=500)
print(TYEX$go_plot)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYEX/cnet_plot.png", width=10, height=10, units="in", res=500)
print(TYEX$cnet_plot)
dev.off()

png("gene_ontology_enrichment_analysis/go_ora_clusterProfiler/TYEX/emap_plot.png", width=10, height=10, units="in", res=500)
print(TYEX$emap_plot)
dev.off()

# =============================================================================
# ## session info
# =============================================================================
sessionInfo()

