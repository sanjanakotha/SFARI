import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import subprocess
from scipy.stats import chisquare
import os

def return_cds_len(folder, code, domain = None):
    loc_df = pd.read_csv("../outputs/mutations/" + folder + "/" + code, sep = "\t", header = None)
    if domain != None:
        loc_df = loc_df[loc_df[3] == domain]
    return sum(loc_df[2] - loc_df[1])

def sort_folder(folder):
    directory = "../outputs/mutations/" + folder 
    files = os.listdir(directory)
    i = 1
    for file in files:  
        if os.path.isdir(directory + "/" + file):
            # skip directories
            continue
        # print(i, file) # gene.bed
        name = file.split(".")[0] # gene
        # print("bedtools sort -i " + directory + "/" + file + " > ../outputs/mutations/cds_bed_format_sorted/" + name + ".bed")
        subprocess.call("bedtools sort -i " + directory + "/" + file + " > " + directory + "/sorted/" + name + ".bed", shell = True)
        i+=1

def generate_row(ENST_code, uniprotID_ENST_mapping_dict, min_var_freq_threshold, max_var_freq_threshold, domain_types, i = 1, print_output = True, variant_filename = "expanded_iWES_v2_variants"):
    tf_lengths, uniprotIDs, AD_lengths, cds_no_syns, AD_no_syns = [],[],[],[],[]
    DBD_lengths, DBD_no_syns = [], []
    lengths = []
    domain_lengths = {}
    domain_counts = {}
    
    for domain_type in domain_types:
        domain_lengths[domain_type] = []
    
    # Building a dataframe
    if print_output:
        print(i, ENST_code)
    i += 1

    # Sort domains
    sort_folder("domains_bed_format")

    # Get uniprotID
    uniprotID = uniprotID_ENST_mapping_dict[ENST_code]
    uniprotIDs.append(uniprotID)

    # Get TF length
    tf_lengths.append(return_cds_len("cds_bed_format", ENST_code))

    # Get sum of domain lengths
    for domain_type in domain_types:
        domain_lengths[domain_type].append(return_cds_len("domains_bed_format", uniprotID, domain_type))

    # Count # rows that are non-syn and below variant threshold in CDS clinvar annotated
    cds_clinvar = pd.read_csv("../outputs/mutations/cds_" + variant_filename + "_snv_classified/" + ENST_code + ".bed", sep = "\t", header = None)
    #display(cds_clinvar)
    cds_no_syn_and_below_count = len(cds_clinvar[(cds_clinvar[cds_clinvar.columns[-1]] == "No-Syn") 
                                                 & (cds_clinvar[cds_clinvar.columns[-4]].astype(float) > min_var_freq_threshold)
                                                 & (cds_clinvar[cds_clinvar.columns[-4]].astype(float) <= max_var_freq_threshold)])
    cds_no_syns.append(cds_no_syn_and_below_count)

    # Count # rows that are non-syn and below variant threshold per domain type
    domains_clinvar = pd.read_csv("../outputs/mutations/domains_" + variant_filename + "_snv_classified/" + ENST_code + ".bed", sep = "\t", header = None)
    #display(domains_clinvar)
    domains_no_syn_and_below = domains_clinvar[(domains_clinvar[domains_clinvar.columns[-1]] == "No-Syn") 
                                               & (domains_clinvar[domains_clinvar.columns[-4]].astype(float) > min_var_freq_threshold)
                                               & (domains_clinvar[domains_clinvar.columns[-4]].astype(float) <= max_var_freq_threshold)]

    for domain_type in domain_types:
        relev_rows = domains_no_syn_and_below[domains_no_syn_and_below[3] == domain_type]
        domain_counts[domain_type] = len(relev_rows)

    results = pd.DataFrame(data = {"uniprotID": uniprotIDs,
                           "TF_cds_length": tf_lengths,
                          "TF_missense" : cds_no_syns
                          })
    results["TF_missense_prop"] = results["TF_missense"] / results["TF_cds_length"]
        
    for domain_type in domain_types:
        results[domain_type + "_cds_length"] = domain_lengths[domain_type]
        results[domain_type + "_missense"] = domain_counts[domain_type]
        results[domain_type + "_missense_prop"] = results[domain_type + "_missense"] / results[domain_type + "_cds_length"]

    return results

def generate_df(ENST_codes, uniprotID_ENST_mapping_dict, min_var_freq_threshold, max_var_freq_threshold, domain_types, print_output = True, variant_filename = "expanded_iWES_v2_variants"):
    i = 1
    rows = []
    for ENST_code in ENST_codes:
        rows.append(generate_row(ENST_code, uniprotID_ENST_mapping_dict, min_var_freq_threshold, max_var_freq_threshold, domain_types, i = i, print_output = print_output, variant_filename = variant_filename))
        i += 1
    return pd.concat(rows).reset_index(drop = True)

import scipy.stats as stats

def add_fisher_p_vals_vs_control(results, domain_type, control = "TF", domain_descrip = ""):
    fisher_exact_p_vals = []

    for i in results.index:
        domain_var_nt = results[domain_type + domain_descrip + "_missense"].iloc[i]
        domain_nt_len = results[domain_type + "_cds_length"].iloc[i]
        
        
        TF_var_nt = results[control + domain_descrip + "_missense"].iloc[i]
        TF_nt_len = results[control + "_cds_length"].iloc[i]

        if control == "TF":
            data = [[domain_var_nt, TF_var_nt - domain_var_nt],
                   [domain_nt_len - domain_var_nt, TF_nt_len - TF_var_nt - domain_nt_len + domain_var_nt]]
        
        else: 
            data = [[domain_var_nt, TF_var_nt],
                    [domain_nt_len - domain_var_nt, TF_nt_len - TF_var_nt]]
            print(data)
        
        # display(results.iloc[i])
        #print(data)

        p_val = stats.fisher_exact(data)[1]
        fisher_exact_p_vals.append(p_val)

    results[domain_type + domain_descrip + "vs" + control + domain_descrip + "_fisher_exact_p_vals"] = fisher_exact_p_vals
    results = results.sort_values(by = domain_type + "vs" + control + "_fisher_exact_p_vals")
    results = results.reset_index(drop = True)
    results = results.reset_index()
    results

def benjamini_hochberg(p_values, alpha):
    """
    Returns decisions on p-values using Benjamini-Hochberg.
    
    Inputs:
        p_values: array of p-values
        alpha: desired FDR (FDR = E[# false positives / # positives])
    
    Returns:
        decisions: binary array of same length as p-values, where `decisions[i]` is 1
        if `p_values[i]` is deemed significant, and 0 otherwise
    """

    # Sorting p-values
    p_vals_df = pd.DataFrame({"p_values" : p_values})
    p_vals_df = p_vals_df.sort_values(by = "p_values")
    p_vals_df = p_vals_df.reset_index()
    p_vals_df = p_vals_df.reset_index()

    # Drawing the line y = k * a / m
    p_vals_df["BH_crit_val"] = alpha * (p_vals_df["level_0"] + 1) / len(p_vals_df)

    # Finding the largest p-value under the line
    if sum(p_vals_df["p_values"] <= p_vals_df['BH_crit_val']) >0:
        threshold = p_vals_df[p_vals_df["p_values"] <= p_vals_df['BH_crit_val']]['p_values'].iloc[-1]
        p_vals_df = p_vals_df.sort_values(by = "index")
        decisions = p_vals_df["p_values"] <= threshold    
        return np.array(decisions)
    else:
        return [False for _ in p_vals_df["p_values"]]

def return_bh_sig(results, col_name, alpha):
    bh_decisions = benjamini_hochberg(results[col_name], alpha)
    return results[bh_decisions]
