# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 09:11:13 2024

@author: Chinaza
"""

import os, pandas as pd, matplotlib.pyplot as plt
from gtfparse import read_gtf
from matplotlib_venn import venn2

print(os.getcwd())

os.chdir("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/")

print(os.getcwd())

# =============================================================================
# import NCBI gtf file and clean it
# =============================================================================

ncbi_refseq_df = read_gtf("large_files/genomic.gtf", result_type='pandas')

ncbi_refseq_filtered_df = ncbi_refseq_df[ncbi_refseq_df['feature'] == 'gene']

ncbi_refseq_filtered_df = ncbi_refseq_filtered_df[['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'gene_id', 'transcript_id', 'db_xref', 'gbkey', 'gene', 'gene_biotype', 'locus_tag']]

ncbi_refseq_filtered_df.reset_index(drop=True, inplace=True)


# =============================================================================
# function links NCBI and Ensembl IDs and compare DEGs from current analysis 
# and nemalife analysis
# =============================================================================

def compare(nemalife,my_analysis,first_figure,second_figure,third_figure):
    
    # =====> NEMALIFE <=====
    nemalife_df = pd.read_csv(nemalife)
    
    # 1st attempt at linking NCBI and Ensembl IDs
    for index, gene_name in nemalife_df["Gene.name"].items():
        for index2, gene_name2 in ncbi_refseq_filtered_df["gene"].items():
            if gene_name == gene_name2:
                nemalife_df.at[index, "ncbi_refseq_gene_id"] = ncbi_refseq_filtered_df.at[index2, 'gene_id']
    
    nan_count = (nemalife_df["ncbi_refseq_gene_id"]== "nan").sum()
    
    if nan_count > 0:    
        print("\nAfter 1st attempt at linking NCBI and Ensembl IDs \n")
        print(f"Genes not yet linked: {nan_count}\n")
        print(nemalife_df[nemalife_df["ncbi_refseq_gene_id"]== "nan"])
        
        # 2nd attempt at linking NCBI and Ensembl IDs
        # Looking at the file closely, some gene name in nemalife_df actually exist in ncbi_refseq_filtered_df but ncbi_refseq_filtered_df has a different gene name; however, the value after ncbi_refseq_filtered_df gene_id "CELE_" contains this name. This could be because ncbi_refseq_filtered_df is a recent annotation (downloaded January 2024) while the annotation source for nemalife_df was an older version. Or this could be because they are from different source (NCBI vs Ensembl).
        for index, gene_name in nemalife_df["Gene.name"].items():
            if nemalife_df.at[index, 'ncbi_refseq_gene_id'] == 'nan':
                for index2, gene_name2 in ncbi_refseq_filtered_df["gene_id"].items():
                    if gene_name in gene_name2:
                        nemalife_df.at[index, "ncbi_refseq_gene_id"] = ncbi_refseq_filtered_df.at[index2, 'gene_id']
        
        nan_count = (nemalife_df["ncbi_refseq_gene_id"]== "nan").sum()
        if nan_count > 0: 
            print("\n\nAfter 2nd attempt at linking NCBI and Ensembl IDs \n")
            print(f"Genes not yet linked: {nan_count}\n")
            print(nemalife_df[nemalife_df["ncbi_refseq_gene_id"]== "nan"])
            
            # 3rd attempt at linking NCBI and Ensembl IDs
            # Linking any leftover genes using "ID" in nemalife and "db_xref" in ncbi_refseq_filtered_df
            for index, ID in nemalife_df["ID"].items():
                if nemalife_df.at[index, 'ncbi_refseq_gene_id'] == 'nan':
                    for index2, db_xref in ncbi_refseq_filtered_df["db_xref"].items():
                        if ID in db_xref:
                            nemalife_df.at[index, "ncbi_refseq_gene_id"] = ncbi_refseq_filtered_df.at[index2, 'gene_id']
                            
            nan_count = (nemalife_df["ncbi_refseq_gene_id"]== "nan").sum()
            if nan_count > 0:               
                print("\n\nAfter 3rd attempt at linking NCBI and Ensembl IDs \n")
                print(f"Genes not yet linked: {nan_count}")                                    
                print(nemalife_df[nemalife_df["ncbi_refseq_gene_id"]== "nan"])
            else:
                print("\nAll genes from Ensembl file linked to NCBI ID after 3rd attempt")  
        else:
            print("\nAll genes from Ensembl file linked to NCBI ID after 2nd attempt")      
    else:
        print("\nAll genes from Ensembl file linked to NCBI ID after 1st attempt")
        
    # =====> MY ANALYSIS <=====
    header = ["ncbi_refseq_gene_id", "baseMean", "log2FoldChange", "lfcSE",	"stat",	"pvalue", "padj"]
    my_analysis = pd.read_table(my_analysis, skiprows=1, names=header)
    
    # =====> COMPARISON AND VENN DIAGRAM <=====
    venn2([set(my_analysis['ncbi_refseq_gene_id']), set(nemalife_df['ncbi_refseq_gene_id'])], set_labels=('current_analysis', 'nemalife_analysis'))
    plt.title('All differentially expressed genes \n(P-adj < 0.05 & absolute(LFC) > 1)')
    plt.savefig(f"data_analysis/data_and_result/compare_DEGs/{first_figure}.png", format="png", dpi=2400, transparent=True)
    plt.show()
   
    DEGs_not_in_nemalife_analysis = set(my_analysis['ncbi_refseq_gene_id']) - set(nemalife_df['ncbi_refseq_gene_id'])
    
    DEGs_not_in_my_analysis = set(nemalife_df['ncbi_refseq_gene_id']) - set(my_analysis['ncbi_refseq_gene_id'])
    
    print("DEGs_not_in_nemalife_analysis:", list(DEGs_not_in_nemalife_analysis))
    print("DEGs_not_in_my_analysis:", list(DEGs_not_in_my_analysis))
    
    
    # =====> UPREGULATED GENES <=====
    nemalife_upregaulated_df = nemalife_df[nemalife_df['log2FoldChange'] > 1]
    my_analysis_upregaulated_df = my_analysis[my_analysis['log2FoldChange'] > 1]
    
    venn2([set(my_analysis_upregaulated_df['ncbi_refseq_gene_id']), set(nemalife_upregaulated_df['ncbi_refseq_gene_id'])], set_labels=('current_analysis', 'nemalife_analysis'))
    plt.title('Upregulated')
    plt.savefig(f"data_analysis/data_and_result/compare_DEGs/{second_figure}.png", format="png", dpi=2400, transparent=True)
    plt.show()
    
    
    # =====> DOWNREGULATED GENES <=====
    nemalife_downregaulated_df = nemalife_df[nemalife_df['log2FoldChange'] < -1]
    my_analysis_downregaulated_df = my_analysis[my_analysis['log2FoldChange'] < -1]
        
    venn2([set(my_analysis_downregaulated_df['ncbi_refseq_gene_id']), set(nemalife_downregaulated_df['ncbi_refseq_gene_id'])], set_labels=('current_analysis', 'nemalife_analysis'))
    plt.title('Downregulated')
    plt.savefig(f"data_analysis/data_and_result/compare_DEGs/{third_figure}.png", format="png", dpi=2400, transparent=True)
    plt.show()


# =============================================================================
# Control vs TPMG
# =============================================================================

nemalife_tpmg = 'some_NemaLife_files/deseq2/Control-vs-TPMG/Significant-DEGs.csv'
my_analysis_tpmg = 'data_analysis/data_and_result/differential_expression_analysis/deseq2/tpmg_up_and_down.txt'

compare(nemalife_tpmg, my_analysis_tpmg, "control_vs_tpmg_all", "control_vs_tpmg_up", "control_vs_tpmg_down")


# =============================================================================
# Control vs TYFD
# =============================================================================

nemalife_tyfd = 'some_NemaLife_files/deseq2/Control-vs-TYFD/Significant-DEGs.csv'
my_analysis_tyfd = 'data_analysis/data_and_result/differential_expression_analysis/deseq2/tyfd_up_and_down.txt'

compare(nemalife_tyfd, my_analysis_tyfd, "control_vs_tyfd_all", "control_vs_tyfd_up", "control_vs_tyfd_down")


# =============================================================================
# Control vs TYEX
# =============================================================================

nemalife_tyex = 'some_NemaLife_files/deseq2/Control-vs-TYEX/Significant-DEGs.csv'
my_analysis_tyex = 'data_analysis/data_and_result/differential_expression_analysis/deseq2/tyex_up_and_down.txt'

compare(nemalife_tyex, my_analysis_tyex, "control_vs_tyex_all", "control_vs_tyex_up", "control_vs_tyex_down")
