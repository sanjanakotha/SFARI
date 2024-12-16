# ---
# jupyter:
#   jupytext:
#     formats: py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
# %load_ext autoreload

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

import AD_predictor_tools
import AD_comparison_tools
import PlottingTools
from Bio.Seq import Seq
# -

uniprotID_ENST_mapping_df = pd.read_csv("../data/SFARI_tf_ENST_codes.csv", index_col = 0)
mart_output = pd.read_csv("../data/mart_export2.txt", sep = "\t")
mart_output = mart_output.dropna(subset = ["Genomic coding start", "Genomic coding end"])


# +
def return_ENST(uniprotID):
    if uniprotID_ENST_mapping_df["uniprotID"].str.contains(uniprotID).any():
        ENST = uniprotID_ENST_mapping_df[uniprotID_ENST_mapping_df["uniprotID"] == uniprotID]["ENST"].iloc[0]
        return ENST
    else:
        print("WARNING: " + uniprotID + " has no corresponding ENST code")
        return ""

# 1-based indexing
def return_mart_output(uniprotID):
    ENST = return_ENST(uniprotID)
    rV = mart_output[mart_output["Transcript stable ID version"] == ENST]
    return  rV.sort_values(by = "Exon rank in transcript")

# Note: by default - 1-based indexing
def return_protein_bed_df(uniprotID, indexing = 1):
    mart_output = return_mart_output(uniprotID)
        
    # Cleaning the mart output
    mart_output["Genomic coding start"] = mart_output["Genomic coding start"].astype(int)
    mart_output["Genomic coding end"] = mart_output["Genomic coding end"].astype(int)
    mart_output = mart_output.rename(columns = {"Chromosome/scaffold name" : "chr",
                                                                "Genomic coding start" : "start",
                                                                "Genomic coding end" : "end"})
    mart_output = mart_output[["chr", "start", "end"]]

    if indexing == 0:
        mart_output["start"] -= 1
        
    return mart_output

# Note: by default - 1-based indexing
def return_domain_bed_df(uniprotID, AA_start, AA_end, indexing = 1):
    mart_output = return_mart_output(uniprotID)
    
    if len(mart_output) == 0:
        return return_protein_bed_df(uniprotID)
    #display(mart_output)
    
    # Converting amino acid start and end to corresponding bp start end
    bp_start = (AA_start - 1) * 3
    bp_end = AA_end * 3
    
    # Cleaning the mart output
    mart_output["Genomic coding start"] = mart_output["Genomic coding start"].astype(int)
    mart_output["Genomic coding end"] = mart_output["Genomic coding end"].astype(int)


    # Only keeping rows that fall within bp start and end range
    filtered_mart_output = mart_output[(mart_output["CDS end"] > bp_start) \
                        & (mart_output["CDS start"] < bp_end)]
    filtered_mart_output = filtered_mart_output.reset_index(drop = True)

    if mart_output["Strand"].iloc[0] == 1:
        
        # Adjust first row
        add = (bp_start - filtered_mart_output["CDS start"].iloc[0] + 1)
        filtered_mart_output.at[0, "Genomic coding start"] += add

        # Adjust last row
        subtract = filtered_mart_output["CDS end"].iloc[-1] - bp_end
        #display(mart_output)
        last_index = filtered_mart_output.index[-1]
        filtered_mart_output.at[last_index, "Genomic coding end"] -= subtract


    elif mart_output["Strand"].iloc[0] == -1:
        
        # Adjust the first row
        subtract = bp_start - filtered_mart_output["CDS start"].iloc[0] + 1
        filtered_mart_output.at[0, "Genomic coding end"] -= subtract

        # Adjust the last row
        add = filtered_mart_output["CDS end"].iloc[-1] - bp_end
        last_index = filtered_mart_output.index[-1]
        filtered_mart_output.at[last_index, "Genomic coding start"] += add

    # Formatting
    filtered_mart_output = filtered_mart_output.sort_values(by = "Genomic coding start")
    filtered_mart_output = filtered_mart_output.rename(columns = {"Chromosome/scaffold name" : "chr",                                                             "Genomic coding start" : "start",
                                                                     "Genomic coding end" : "end"})
    filtered_mart_output = filtered_mart_output[["chr", "start", "end"]]

    if indexing == 0:
        filtered_mart_output["start"] -= 1
        
    return filtered_mart_output

def translate(seq):
    seq_obj = Seq(seq)
    return str(seq_obj.translate())

# Note: by default - 1-based indexing
def AAseq_using_genomic_coords(uniprotID, AA_start, AA_end):
    mart_output = return_mart_output(uniprotID)
    bed_df = return_domain_bed_df(uniprotID, AA_start, AA_end)
    #display(bed_df)
    
    if len(bed_df) == 0:
        print("WARNING: " + uniprotID + " has no corresponding ENST code, so can't get AA sequence.")
        return ""
        
    # Source for following lines: https://stackoverflow.com/questions/58010068/how-to-explode-range-from-two-column-to-rows
    zipped = zip(bed_df['chr'],
                 bed_df['start'], 
                bed_df['end'])
    genomic_coord_df = pd.DataFrame([(i, y) for i,j,k in zipped for y in range(j, k+1)], columns=['chr', 'position'])
    indices = list(genomic_coord_df["position"] - 1)
    
    chromosome = bed_df["chr"].iloc[0]
    chr_df = AD_predictor_tools.makeFullLengthProteinDF("../data/hg38/chr" + chromosome + ".fa");

    DNA_seq = ""
    for index in indices:
        DNA_seq += chr_df["AAseq"].iloc[0][index]
        
    # print("DNA seq:")
    # print(DNA_seq)
    # print()
    
    if mart_output["Strand"].iloc[0] == 1:
        return translate(DNA_seq)
    else:
        seq_obj = Seq(DNA_seq)
        return str(seq_obj.reverse_complement().translate())
    
def return_domain_bed_df_for_full_prot(uniprotID, indexing = 1, known_ADs_filepath = "../output/SFARI_TF_known_ADs.csv"):
    known_ADs = pd.read_csv("../output/SFARI_TF_known_ADs.csv", index_col = 0)
    
    ADs = known_ADs[known_ADs["uniprotID"] == uniprotID]
    coord_dfs = []
    
    for start, end in zip(ADs["Start"], ADs["End"]):
        coord_dfs.append(return_domain_bed_df(uniprotID, start, end, indexing = 0))
    
    return pd.concat(coord_dfs)
# +
def save_bed_from_df(df, filepath, add_chr = True):
    if add_chr:
        df["chr"] = "chr" + df["chr"]
    df.to_csv(filepath, index=False, sep='\t', header=None)
    
def save_protein_bed(uniprotID, filepath):
    df = return_protein_bed_df(uniprotID, indexing = 0)
    save_bed_from_df(df, filepath)

def save_domain_bed(uniprotID, start, end, filepath):
    df = return_domain_bed_df(uniprotID, start, end, indexing = 0)
    save_bed_from_df(df, filepath)
    
def save_domains_bed(uniprotID, filepath, known_ADs_filepath = "../output/SFARI_TF_known_ADs.csv"):
    df = return_domain_bed_df_for_full_prot(uniprotID, indexing = 0, known_ADs_filepath = known_ADs_filepath)
    save_bed_from_df(df, filepath)
# -


