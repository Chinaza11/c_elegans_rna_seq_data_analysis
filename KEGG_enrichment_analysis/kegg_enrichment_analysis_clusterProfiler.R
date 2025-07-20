# =============================================================================
# ## load packages
# =============================================================================

library(clusterProfiler)
library(tibble)
library(pathview)
library(ggplot2)

# =============================================================================
# ## function
# =============================================================================
set.seed(1234)

KEGG_ORA_n_GSEA_analysis = function(path,figures_path_ora,figures_path_gsea){

        setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/")
        
        modify_plot = theme(plot.title = element_text(size = 24, face ="bold", hjust = 0.5))
        
        df = read.table(path, header=TRUE, sep="\t")
        df = rownames_to_column(df, var = "ncbi_locus_tag")
        #print(head(df))
        
        ids = read.csv("DAVID-locus_tag_to_entrez_id.csv")
        names(ids)[names(ids) == "From"] = "ncbi_locus_tag"
        names(ids)[names(ids) == "To"] = "Entrez_id"
        #print(head(ids))
        
        merged_df = merge(df, ids, by = "ncbi_locus_tag", all = TRUE)
        #print(head(merged_df))
        
        genes_of_interest = as.character((subset(merged_df, padj < 0.05 & abs(log2FoldChange) > 1))$ncbi_locus_tag)
        
                                          
        # KEGG pathway: over-representation analysis using Fisher's test
        kk_ora = enrichKEGG(gene         = genes_of_interest,
                        organism     = 'cel',
                        universe      = merged_df$ncbi_locus_tag,
                        pvalueCutoff = 0.05)
        #head(kk_ora)
        result_kk_ora = data.frame(kk_ora)
        
        dotplot_kk_ora = if((dim(result_kk_ora)[1]) > 0){
          dotplot(kk_ora, showCategory = dim(result_kk_ora)[1]) + 
            labs(title="Enriched Pathway", 
                 caption="GeneRatio ==> ratio of differential expressed genes that are annotated in a pathway\nCount ==> GeneCount") + 
            modify_plot
        }
        
        # prepare geneList for GSEA
        geneList = merged_df$log2FoldChange
        names(geneList) = merged_df$ncbi_locus_tag
        geneList = sort(geneList, decreasing = TRUE)
        
        
        # KEGG pathway: GSEA using KS test
        kk_gsea = gseKEGG(geneList     = geneList,
                       organism     = 'cel',
                       pvalueCutoff = 0.05,
                       verbose      = FALSE,
                       seed = TRUE)
        #head(kk_gsea)
        result_kk_gsea = data.frame(kk_gsea)
        
        dotplot_kk_gsea =  if((dim(result_kk_gsea)[1]) > 0){
          dotplot(kk_gsea, showCategory = dim(result_kk_gsea)[1]) + 
            labs(title="Enriched Pathway", 
                 caption="GeneRatio ==> ratio of differential expressed genes that are annotated in a pathway\nCount ==> GeneCount") + 
            modify_plot
        }
        
        # Visualize enriched KEGG pathways
        get_pathview = function(pathway_id, gene_data, out, figures_path){
          setwd(figures_path)  
          pathview(gene.data  = gene_data,
                            pathway.id = pathway_id,
                            gene.idtype = "KEGG",
                            species    = "cel",
                            out.suffix = out)}
        
        # GSEA
        for (i in result_kk_gsea$ID){
          print(i)
          get_pathview(i,geneList,"GSEA",figures_path_gsea)}
        
        # ORA
        for (i in result_kk_ora$ID){
          print(i)
          get_pathview(i,geneList,"ORA",figures_path_ora)}
        
        return(output=list(result_kk_ora = result_kk_ora, 
                           dotplot_kk_ora = dotplot_kk_ora,
                           result_kk_gsea = result_kk_gsea,
                           dotplot_kk_gsea = dotplot_kk_gsea))
}

# =============================================================================
# ## TPMG vs Control
# =============================================================================

path = "differential_expression_analysis/deseq2/tpmg_vs_control_df.txt"
figures_path_ora="G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TPMG/ORA-Fisher_test"
figures_path_gsea="G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TPMG/GSEA-KS_test"

TPMG = KEGG_ORA_n_GSEA_analysis(path, figures_path_ora, figures_path_gsea)


setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TPMG")

#write.table(TPMG$result_kk_ora, "ORA-Fisher_test/result_kk_ora.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("ORA-Fisher_test/dotplot_kk_ora.png", width=12, height=12, units="in", res=500)
print(TPMG$dotplot_kk_ora)
dev.off()

#write.table(TPMG$result_kk_gsea, "GSEA-KS_test/result_kk_gsea.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("GSEA-KS_test/dotplot_kk_gsea.png", width=12, height=12, units="in", res=500)
print(TPMG$dotplot_kk_gsea)
dev.off()

# =============================================================================
# ## TYFD vs Control
# =============================================================================

path = "differential_expression_analysis/deseq2/tyfd_vs_control_df.txt"
figures_path_ora="G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TYFD/ORA-Fisher_test"
figures_path_gsea="G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TYFD/GSEA-KS_test"

TYFD = KEGG_ORA_n_GSEA_analysis(path, figures_path_ora, figures_path_gsea)


setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TYFD")

#write.table(TYFD$result_kk_gsea, "GSEA-KS_test/result_kk_gsea.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("GSEA-KS_test/dotplot_kk_gsea.png", width=12, height=12, units="in", res=500)
print(TYFD$dotplot_kk_gsea)
dev.off()

# =============================================================================
# ## TYEX vs Control
# =============================================================================

path = "differential_expression_analysis/deseq2/tyex_vs_control_df.txt"
figures_path_ora="G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TYEX/ORA-Fisher_test"
figures_path_gsea="G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TYEX/GSEA-KS_test"

TYEX = KEGG_ORA_n_GSEA_analysis(path, figures_path_ora, figures_path_gsea)


setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TYEX")

#write.table(TYEX$result_kk_ora, "ORA-Fisher_test/result_kk_ora.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("ORA-Fisher_test/dotplot_kk_ora.png", width=12, height=12, units="in", res=500)
print(TYEX$dotplot_kk_ora)
dev.off()

#write.table(TYEX$result_kk_gsea, "GSEA-KS_test/result_kk_gsea.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

png("GSEA-KS_test/dotplot_kk_gsea.png", width=12, height=12, units="in", res=500)
print(TYEX$dotplot_kk_gsea)
dev.off()

# =============================================================================
# ## session info
# =============================================================================
sessionInfo()