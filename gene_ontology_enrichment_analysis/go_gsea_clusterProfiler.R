# =============================================================================
# ## setting working directory and loading packages
# =============================================================================
setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result")

library(clusterProfiler)
library(org.Ce.eg.db)
library(tibble)
library(enrichplot)
library(ggplot2)


set.seed(1234)


# =============================================================================
# ## function
# =============================================================================

gsea_analysis = function(path,treatment){

    df = read.table(path, header=TRUE, sep="\t")
    df = rownames_to_column(df, var = "ncbi_locus_tag")
    print(head(df))
    
    ids = read.csv("DAVID-locus_tag_to_entrez_id.csv")
    names(ids)[names(ids) == "From"] = "ncbi_locus_tag"
    names(ids)[names(ids) == "To"] = "Entrez_id"
    print(head(ids))
    
    merged_df = merge(df, ids, by = "ncbi_locus_tag", all = TRUE)
    print(head(merged_df))
    
    
    # prepare geneList for GSEA
    geneList = merged_df$log2FoldChange
    names(geneList) = merged_df$Entrez_id
    geneList = sort(geneList, decreasing = TRUE)
    
    #TPMG ==> nPermSimple = 100000, eps = 0
    #TYFD ==> no modifications
    #TYEX ==> nPermSimple = 10000, eps = 0
    
    if (treatment == "TPMG"){
      eps = 0
      nPermSimple = 100000
    } else if (treatment == "TYFD") {
      eps = 1e-10
      nPermSimple = 1000
    } else if (treatment == "TYEX") {
      eps = 0
      nPermSimple = 10000
      }
    
    ego <- gseGO(geneList     = geneList,
                 OrgDb        = org.Ce.eg.db,
                 ont          = "BP",
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 seed         = TRUE,
                 eps          = eps,
                 nPermSimple  = nPermSimple)
    result = data.frame(ego)
    
    
    ## use simplify to remove redundant terms
    ego2 = simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
    result2 = data.frame(ego2)
    
    return(list(df_all_terms = result, df_redudant_terms_removed = result2))

}

# =============================================================================
# ## Control vs TPMG 
# =============================================================================

path = "differential_expression_analysis/deseq2/tpmg_vs_control_df.txt"
treatment = "TPMG"

tpmg_result = gsea_analysis(path,treatment)

view(tpmg_result$df_all_terms)
view(tpmg_result$df_redudant_terms_removed)

#write.table(tpmg_result$df_all_terms, "gene_ontology_enrichment_analysis/go_gsea_clusterProfiler/TPMG/result_all_terms.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(tpmg_result$df_redudant_terms_removed, "gene_ontology_enrichment_analysis/go_gsea_clusterProfiler/TPMG/result_redundant_terms_removed.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

# =============================================================================
# ## Control vs TYFD
# =============================================================================

path = "differential_expression_analysis/deseq2/tyfd_vs_control_df.txt"
treatment = "TYFD"

tyfd_result = gsea_analysis(path,treatment)

view(tyfd_result$df_all_terms)
view(tyfd_result$df_redudant_terms_removed)

#write.table(tyfd_result$df_all_terms, "gene_ontology_enrichment_analysis/go_gsea_clusterProfiler/TYFD/result_all_terms.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(tyfd_result$df_redudant_terms_removed, "gene_ontology_enrichment_analysis/go_gsea_clusterProfiler/TYFD/result_redundant_terms_removed.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

# =============================================================================
# ## Control vs TYEX
# =============================================================================

path = "differential_expression_analysis/deseq2/tyex_vs_control_df.txt"
treatment = "TYEX"

tyex_result = gsea_analysis(path,treatment)

view(tyex_result$df_all_terms)
view(tyex_result$df_redudant_terms_removed)

#write.table(tyex_result$df_all_terms, "gene_ontology_enrichment_analysis/go_gsea_clusterProfiler/TYEX/result_all_terms.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(tyex_result$df_redudant_terms_removed, "gene_ontology_enrichment_analysis/go_gsea_clusterProfiler/TYEX/result_redundant_terms_removed.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")
# =============================================================================
# ## session info
# =============================================================================

sessionInfo()

