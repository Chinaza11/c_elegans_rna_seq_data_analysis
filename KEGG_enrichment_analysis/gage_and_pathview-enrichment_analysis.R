# =============================================================================
# ## setting working directory and loading packages
# =============================================================================

library(gage)
library(pathview)
library(ggplot2)

# =============================================================================
# ## function
# =============================================================================

gage_analysis = function(path, p_value, title, subtitle, pathview_path){
  home = "G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/"
  setwd(home)
  
    #========> files prep <========
    df = read.table(path, header=TRUE, sep="\t")
    df = rownames_to_column(df, var = "ncbi_locus_tag")
    #print(head(df))
    
    ids = read.csv("DAVID-locus_tag_to_entrez_id.csv")
    names(ids)[names(ids) == "From"] = "ncbi_locus_tag"
    names(ids)[names(ids) == "To"] = "Entrez_id"
    #print(head(ids))
    
    merged_df = merge(df, ids, by = "ncbi_locus_tag", all = TRUE)
    #print(head(merged_df))
    
    exp.fc = merged_df$log2FoldChange
    names(exp.fc) = merged_df$ncbi_locus_tag
    
    kg.cel = kegg.gsets("cel")
    kegg.gs = kg.cel$kg.sets[kg.cel$sigmet.idx]
    #save(kegg.gs, file="KEGG_enrichment_analysis/kegg.cel.sigmet.gsets.RData")
    #load(file = "KEGG_enrichment_analysis/kegg.cel.sigmet.gsets.RData")
    
    #========> GAGE <========
    fc.kegg.p = gage(exp.fc, gsets=kegg.gs, ref = NULL, samp = NULL)
    sel = fc.kegg.p$greater[, "q.val"] < 0.05 & !is.na(fc.kegg.p$greater[, "q.val"])
    path.ids = rownames(fc.kegg.p$greater)[sel]
    sel.l = fc.kegg.p$less[, "q.val"] < 0.05 & !is.na(fc.kegg.p$less[,"q.val"])
    path.ids.l = rownames(fc.kegg.p$less)[sel.l]
    path.ids2 = substr(c(path.ids, path.ids.l), 1, 8)
    print(paste(length(path.ids2), "significant pathways"))
    
    if (length(path.ids2) > 0) {
      greater = data.frame(fc.kegg.p$greater)
      greater_sub = greater[greater$q.val < p_value & !is.na(greater$q.val), ]
      
      less = data.frame(fc.kegg.p$less)
      less_sub = less[less$q.val < p_value & !is.na(less$q.val), ]
      
      combined_df = rbind(greater_sub,less_sub)
    
      barplot = ggplot(combined_df, aes(x=rownames(combined_df), y= -log(q.val))) + 
                  geom_bar(stat="identity",fill="steelblue") + 
                  labs(title=title, subtitle=subtitle) + 
                  xlab("Pathways") + 
                  scale_y_continuous(limits=c(0,max(-log(combined_df$q.val))+2)) + 
                  coord_flip() +
                  theme_minimal() + 
                  theme(plot.title = element_text(size = 26, face ="bold"),
                        plot.subtitle = element_text(size = 14),
                        axis.title.y = element_text(size=14, face="bold"),
                        axis.title.x = element_text(size=14, face="bold"),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 12))
      
      #========> PATHVIEW <========
      exp.fc2 = merged_df$log2FoldChange
      names(exp.fc2) = merged_df$Entrez_id
      setwd(pathview_path)
      pv.out.list = sapply(TPMG$path.ids2, function(pid) pathview(gene.data = exp.fc2, pathway.id = pid, species = "cel", out.suffix="GAGE"))
      setwd(home)
  }
      
    if (length(path.ids2) > 0) {
      return(list(barplot = barplot,
                  path.ids2 = path.ids2,
                  fc.kegg.p_greater = fc.kegg.p$greater,
                  fc.kegg.p_less = fc.kegg.p$less))
    } else {
      return(list(fc.kegg.p_greater = fc.kegg.p$greater,
                  fc.kegg.p_less = fc.kegg.p$less))
    }
}

# =============================================================================
# ## TPMG vs Control
# =============================================================================

path = "differential_expression_analysis/deseq2/tpmg_vs_control_df.txt"
p_value = 0.05
title = "Control vs TPMG (GAGE)"
subtitle = "Pathways with q-value less than 0.05"
pathview_path = "G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TPMG/GAGE"

TPMG = gage_analysis(path, p_value, title, subtitle, pathview_path)

png("KEGG_enrichment_analysis/TPMG/GAGE/Control_vs_TPMG_0.05.png", 
    width=12, height=12, units="in", res=500)
print(TPMG$barplot)
dev.off()

#write.table(TPMG$fc.kegg.p_greater, "KEGG_enrichment_analysis/TPMG/GAGE/fc.kegg.p_greater.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(TPMG$fc.kegg.p_less, "KEGG_enrichment_analysis/TPMG/GAGE/fc.kegg.p_less.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

# =============================================================================
# ## TYFD vs Control
# =============================================================================

path = "differential_expression_analysis/deseq2/tyfd_vs_control_df.txt"
p_value = 0.05
title = "Control vs TYFD (GAGE)"
subtitle = "Pathways with q-value less than 0.05"
pathview_path = "G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TYFD/GAGE"

TYFD = gage_analysis(path, p_value, title, subtitle, pathview_path)

# no enriched pathway

#write.table(TYFD$fc.kegg.p_greater, "KEGG_enrichment_analysis/TYFD/GAGE/fc.kegg.p_greater.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(TYFD$fc.kegg.p_less, "KEGG_enrichment_analysis/TYFD/GAGE/fc.kegg.p_less.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

# =============================================================================
# ## TYEX vs Control
# =============================================================================

path = "differential_expression_analysis/deseq2/tyex_vs_control_df.txt"
p_value = 0.05
title = "Control vs TYEX (GAGE)"
subtitle = "Pathways with q-value less than 0.05"
pathview_path = "G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/KEGG_enrichment_analysis/TYEX/GAGE"

TYEX = gage_analysis(path, p_value, title, subtitle, pathview_path)

# no enriched pathway

#write.table(TYEX$fc.kegg.p_greater, "KEGG_enrichment_analysis/TYEX/GAGE/fc.kegg.p_greater.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

#write.table(TYEX$fc.kegg.p_less, "KEGG_enrichment_analysis/TYEX/GAGE/fc.kegg.p_less.txt", quote=FALSE, col.names=TRUE, row.names=TRUE, sep="\t")

# =============================================================================
# ## session info
# =============================================================================

sessionInfo()