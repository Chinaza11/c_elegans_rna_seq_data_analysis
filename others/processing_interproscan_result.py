# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 15:20:39 2024

@author: Chinaza
"""

import os, pandas as pd

print(os.getcwd())

os.chdir("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/large_files/")

print(os.getcwd())

# =============================================================================
# 
# =============================================================================

header = ["Protein accession", "Sequence MD5 digest", "Sequence length", "Analysis", "Signature accession", "Signature description", "Start location", "Stop location", "Score or e-value", "Status", "Date", "InterPro annotations - accession", "InterPro annotations - description", "GO annotations with their source(s)", "Pathways annotations"]

interproscan_result = pd.read_table("protein.faa.tsv", names=header)
interproscan_result.head()

unique_protein_accession  = interproscan_result['Protein accession'].unique()

# =============================================================================
# 
# =============================================================================

annot_df = pd.DataFrame()

for index, go_term in interproscan_result["GO annotations with their source(s)"].items():
    if pd.isna(go_term) != True and interproscan_result.at[index, "GO annotations with their source(s)"] != "-": 
        print(go_term)
        annot_df.at[index, "Protein accession"] = interproscan_result.at[index, "Protein accession"]
        annot_df.at[index, "GO term"] = interproscan_result.at[index, "GO annotations with their source(s)"]

annot_df_1 = annot_df.drop_duplicates()
annot_df_2 = annot_df_1.groupby('Protein accession').agg(lambda x: ', '.join(set(x)))
annot_df_2 = annot_df_2.reset_index()

#14181 proteins are annotated with GO term

