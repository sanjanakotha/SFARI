import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import seaborn as sns
import numpy as np
import re
import metapredict as meta

import warnings
warnings.filterwarnings("ignore")
from statsmodels.nonparametric.smoothers_lowess import lowess

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

import protfasta
import glob
import itertools

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Helvetica" #somethings this one doesnt work
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['figure.dpi'] = 300

# Define colors for the domains
AD_color = 'orange'
RD_color = 'red'
DBD_color = 'purple'

TFs_tbl = pd.read_csv("../soto_analysis/outputs/all_TFs_table_proteins.txt", sep="\t", index_col=0)
lambert_TFs = pd.read_csv("../output/lambert_TFs_10-21-24_with_DBD_coords.csv")
lambert_TFs["uniprotID"] = lambert_TFs["id"].str.split("|").str[1]
lambert_TFs["len"] = lambert_TFs["ProteinSeq"].str.len()
TFs_tbl = pd.merge(TFs_tbl, lambert_TFs, on = 'uniprotID')

def return_new_ranges(col):
    split = TFs_tbl[col].astype(str).str.split(',')
    new_ranges = []
    
    for TF_ranges in split:
        # Handle missing or invalid entries
        if TF_ranges is None or TF_ranges[0] in ["nan", "NA", "None", ""]:
            new_ranges.append([])
            continue
        
        new_TF_ranges = []
        for indiv_range in TF_ranges:
            try:
                start, end = indiv_range.split("-")
                new_TF_ranges.append((int(start), int(end)))
            except (ValueError, AttributeError):
                # Skip bad entries like "NA" or malformed strings
                continue
                
        new_ranges.append(new_TF_ranges)
    
    return new_ranges

TFs_tbl["DBD_ranges"] = return_new_ranges("DBD_coords")
TFs_tbl["AD_ranges"] = return_new_ranges("AD_coords")
TFs_tbl["RD_ranges"] = return_new_ranges("RD_coords")

TFs_tbl[TFs_tbl["uniprotID"].str.contains("Q13950", na=False)]

import itertools

def return_intervals(a):
    a = sorted(set(a))
    b = []
    for k, g in itertools.groupby(enumerate(a), 
        key=lambda t: t[1] - t[0]):
        g = list(g)
        b.append([g[0][1], g[-1][1]]) 
    
    return b
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

def threshold_plot(data, x, y, legend, ax, lw, color_below, color_above, threshold, linestyle = "solid"):
    if data is None:
        data = pd.DataFrame({'x' : x, 
                             'y' : y})
        x = 'x'
        y = 'y'

    if color_below == "alpha":
        sns.lineplot(data = data, x = x, y = y, legend = legend, ax = ax, lw = lw, color = color_above, alpha = 0.4, linestyle = linestyle)
    else:
        sns.lineplot(data = data, x = x, y = y, legend = legend, ax = ax, lw = lw, color = color_below, linestyle = linestyle)

    data = data.reset_index(drop = True)
    
    above = data[(data[y] >= threshold) | pd.isnull(data[y])]
    #display(above)
    
    intervals_above = return_intervals(list(above.index))
    #print(intervals_above)
    #display(data.loc[8])
    for interval in intervals_above:
        #print(interval)
        interval_df = data.loc[interval[0] : interval[1]]
        #display(interval_df)
        sns.lineplot(data = interval_df, x = x, y = y, legend = legend, ax = ax, lw = 1, color = color_above, linestyle = linestyle)

def add_custom_legend(color_dict, ax, bbox_to_anchor = (1, 1)):
    custom_lines = [Line2D([0], [0], markersize=2, color=c, lw=4) for c in color_dict.values()]
    ax.legend(custom_lines, color_dict.keys(), bbox_to_anchor = bbox_to_anchor, frameon = False, fontsize = 'x-small', handlelength = 0.5)

def add_custom_legend_arrows(top_color, top_label, bottom_color, bottom_label, ax, center = 0.5, length = 0.2, padding = 0.05):
    #custom_lines = [Line2D([0], [0], markersize=2, color=c, lw=4) for c in color_dict.values()]

        # Create custom arrows for the legend
    custom_lines = []
    
    # Top arrow (points up)
    ax.add_patch(FancyArrowPatch(
            (1.03, center + padding), 
            (1.03, center + padding + length),  # This makes the arrow point up
            mutation_scale=10, 
            color=top_color, 
            linewidth=1,
            transform=ax.transAxes, clip_on = False
        ))
    
    ax.text(x = 1.055, y = 0.5 * (2 * center + 2 * padding + length), s= top_label, transform=ax.transAxes, fontsize = "x-small", va = "center")
             
    # Bottom arrow (points down)
    ax.add_patch(FancyArrowPatch(
            (1.03, center - padding), 
            (1.03, center - padding - length),  # This makes the arrow point down
            mutation_scale=10, 
            color=bottom_color, 
            linewidth=1,
            transform=ax.transAxes, clip_on = False
        ))
    
    ax.text(x = 1.055, y = 0.5 * (2 * center - 2 * padding - length), s= bottom_label, transform=ax.transAxes, fontsize = "x-small", va = "center")

def plot_rectangle(start, end, color, box_height, ax, y = 0, alpha = 0.7):
    rectangle = patches.Rectangle((start, y), width=end - start, height=box_height, linewidth=1, edgecolor='none', facecolor=color, alpha=alpha)
    ax.add_patch(rectangle)

def merge_intervals(intervals, threshold=10):
    # Sort intervals by their starting value
    intervals.sort(key=lambda x: x[0])
    
    merged = []
    for interval in intervals:
        if not merged:
            merged.append(interval)
        else:
            # Check if the current interval can be merged with the last one in the merged list
            if interval[0] - merged[-1][1] < threshold:
                # Merge the intervals by extending the end of the last interval
                merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
            else:
                # Otherwise, just add the current interval as a new entry
                merged.append(interval)
    
    return merged
    
def plot_domain_rectangles(DBD_coords, AD_coords, RD_coords, box_height, DBD_y, AD_y, RD_y, ax, RD = True, AD_text = True, alpha = 0.7, other_domains = None):

    DBD_coords = merge_intervals(DBD_coords)
    
    if other_domains is not None:
    
        for name, start, end in zip(other_domains["domain_description"], other_domains["start"], other_domains["end"]):
            plot_rectangle(start, end, "gray", box_height, ax, y = 0, alpha = alpha)
        
            if AD_text:
                ax.text(np.mean([start, end]), 0.5 * box_height, 
                         name, 
                         fontsize='small', 
                         color="black", 
                         horizontalalignment='center',
                        verticalalignment='center')

                #gray_track_y_adj = 0
                

        # AD_y -= gray_track_y_adj
        # DBD_y -= gray_track_y_adj
        # RD_y -= gray_track_y_adj
    
    for start, end in AD_coords:
        #print(start, end)
        plot_rectangle(start, end, AD_color, box_height, ax, y = AD_y, alpha = alpha)

        if AD_text:
            # cc_names_row = cc_names[(cc_names["end"] == end) & (cc_names["start"] == start)]

            # if len(cc_names_row) > 0:
            #     AD_name = cc_names_row["Gene Name"].iloc[0]
            #     AD_name = AD_name.split("_")[-1]
            #     # if AD_name[-1] in ["1", "2", "3"]:
            #     #     AD_name = "AD" + AD_name[-1]
            ax.text(np.mean([start, end]), AD_y + 0.5 * box_height, 
                        "AD", 
                        fontsize='small', 
                        color="black", 
                        horizontalalignment='center',
                    verticalalignment='center')

    if RD:
        for start, end in RD_coords:
            if (start != 0) & (end != 0):
            #print(start, end)
            
                plot_rectangle(start, end, RD_color, box_height, ax, y = RD_y , alpha = alpha)
    
                if AD_text:
                    ax.text(np.mean([start, end]), RD_y + 0.5 * box_height, 
                                 "RD", 
                                 fontsize='small', 
                                 color="black", 
                                 horizontalalignment='center',
                                verticalalignment='center')
        
    for start, end in DBD_coords:
        plot_rectangle(start, end, DBD_color, box_height, ax, y = DBD_y, alpha = alpha)

        if AD_text:
            ax.text(np.mean([start, end]), DBD_y + 0.5 * box_height, 
                     "DBD", 
                     fontsize='small', 
                     color="black", 
                     horizontalalignment='center',
                    verticalalignment='center')
        

sns.set_context('paper')

import matplotlib.patches as patches

# Plots DBD, AD, and RD annotations
def plot_annots(uniprotID, ax, box_height = 0.3, RD = False, other_domains = None):
    row = TFs_tbl[TFs_tbl["uniprotID"].str.contains(uniprotID)].iloc[0]

    TF_len = len(row['ProteinSeq'])
    #ENST = row["ENST"]
    DBD_coords = row["DBD_ranges"]
    AD_coords = row["AD_ranges"]
    RD_coords = row["RD_ranges"]

    #plt.figure(figsize = (8, 0.75))
    #ax.axhline(y = 0, color = 'black', linestyle = '-', lw = 2) 
    #ax.axhline(y = box_height, color = 'black', linestyle = '--', lw = 1, alpha = 0.1) 
    #ax.axhline(y = 2 * box_height, color = 'black', linestyle = '--', lw = 1, alpha = 0.1) 
    
    
    plot_domain_rectangles(DBD_coords, 
                           AD_coords, 
                           RD_coords, 
                           box_height,
                          DBD_y = 0.5,
                            AD_y = 0.75,
                          RD_y = 0.25, ax =ax, alpha = 0.7, RD = RD, other_domains = other_domains, AD_text=True)
    
        
    # Adjust plot limits so that the rectangles are visible
    ax.set_ylim(0, 1)  # y-limits remain from 0 to 1, as the plot height is small
    #plt.axis('off')
    
    #plt.yticks([box_height / 2, box_height * 1.5], ['DBD', 'AD'])
    ax.set_yticks([])
    
    y_labels = ax.get_yticklabels()
    
    color_dict = {"AD" : AD_color, "RD" : RD_color, "DBD" : DBD_color}

    # for i, label in enumerate(y_labels):
    #     domain_type = re.search("'(.*)'", str(label))[0][1:-1]
    #     label.set_color(color_dict[domain_type])

    ax.set_xlim(0, TF_len)  # Adjust x-limits to fit rectangles (based on start and end values)
    ax.set_ylim(0, 0.75 + box_height)
    ax.set_xlabel("Protein Residue", labelpad = 0)
    
    # ax2 = ax.twiny()
    # ax2.xaxis.set_ticks_position("bottom")
    # ax2.xaxis.set_label_position("bottom")
    # new_tick_locations = np.arange(0, TF_len * 3 + 1, 300)
    # ax2.spines["bottom"].set_position(("axes", -1))
    
    # ax2.spines["bottom"].set_color(sns.color_palette('tab10')[2])
    
    # #ax2.set_frame_on(True)
    # #ax2.patch.set_visible(False)
    
    # # for sp in ax2.spines.itervalues():
    # #     sp.set_visible(False)
    # #ax2.spines["bottom"].set_visible(True)
    
    # ax2.set_xticks(new_tick_locations)
    # #ax2.set_xticklabels(new_tick_locations)
    # ax2.tick_params(axis='x', colors=sns.color_palette('tab10')[2])
    # ax2.set_xlabel(r"CDS Nucleotide", color = sns.color_palette('tab10')[2], labelpad = 0)
    # sns.despine(left = True, ax = ax2)

    plot_domain_rectangles(DBD_coords, 
                           AD_coords,
                           RD_coords, 1, 
                           RD = False, AD_text = False, alpha = 0.15,
                          DBD_y = 0., 
                           AD_y = 0, 
                           RD_y = 0, 
                           ax = ax)

    sns.despine(left = True, ax = ax)

def plot_disorder(uniprotID, ax):
    row = TFs_tbl[TFs_tbl["uniprotID"].str.contains(uniprotID)].iloc[0]
    #ENST = row["ENST"]
    DBD_coords = row["DBD_ranges"]
    AD_coords = row["AD_ranges"]
    RD_coords = row["RD_ranges"]

    disorder = meta.predict_disorder(row["ProteinSeq"])
    ax.axhline(0.5, color = "gray", linestyle = "--", alpha = 0.7)

    #sns.lineplot(x = np.arange(len(disorder)), y = disorder, ax = ax, color = sns.color_palette('colorblind')[7])

    data = pd.DataFrame({'pos' : np.arange(len(disorder)), 'disorder': disorder})
    
    threshold_plot(data = data, 
                   x = 'pos', 
                   y = 'disorder', 
                   legend = True, 
                   ax = ax, 
                   lw = 1, 
                   color_below = "silver", 
                   color_above = sns.color_palette('colorblind')[2], 
                   threshold = 0.5)
    
    ax.set_xlim(0, len(disorder))


    plot_domain_rectangles(DBD_coords, 
                    AD_coords, 
                    RD_coords, 1, RD = False, AD_text = False, alpha = 0.15,
                              DBD_y = 0, AD_y = 0, RD_y = 0, ax = ax)

    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Disorder")
    #plt.title(cc_TFs["Gene"].iloc[i])
    sns.despine()



    # ax.text(max(data['pos']), 0.52, 'Disordered', ha = "left", va = "bottom", color = sns.color_palette('colorblind')[2])
    # ax.text(max(data['pos']), 0.48, 'Folded', ha = "left", va = "top", color = "gray")

    #plt.show()

def plot_cds_trace(uniprotID, ax, x_axis_spacing = 600, alignment = 'zoonomia'):
    cds_phylo_P = []

    cds_paths = glob.glob(f"../soto_analysis/outputs/mutations/cds_{alignment}_all_TF_cds/*")
    for path in cds_paths:
        ENST = path.split("/")[-1].split(".bed")[0]
        df = pd.read_csv(path, sep = "\t", header = None)
        cds_phylo_P.append(df)

    cds_phylo_P = pd.concat(cds_phylo_P)
    cds_phylo_P = cds_phylo_P.rename(columns = {3: "ENST", 8: "PhyloP"})
    cds_phylo_P

    row = TFs_tbl[TFs_tbl["uniprotID"].str.contains(uniprotID)].iloc[0]
    ENST = row["ENST"]
    DBD_coords = row["DBD_ranges"]
    AD_coords = row["AD_ranges"]
    RD_coords = row["RD_ranges"]

    ENST_phylo_P = cds_phylo_P[cds_phylo_P["ENST"] == ENST]
    
    if ENST_phylo_P[4].iloc[0] == -1:
        ascending = False
    else: 
        #print("pos")
        ascending = True
        
    ENST_phylo_P = ENST_phylo_P.sort_values(by = 1, ascending = ascending)
    ENST_phylo_P["cds_nt"] = np.arange(len(ENST_phylo_P)) / 3

    smoothed = lowess(ENST_phylo_P['PhyloP'], 
                      ENST_phylo_P["cds_nt"], 
                      frac=0.01)  # Adjust frac for smoothness
    smoothed_df = pd.DataFrame(smoothed, columns=["cds_nt", 'smoothed'])

    # sns.scatterplot(data=ENST_phylo_P, 
    #                 x="cds_nt", 
    #                 y='PhyloP', 
    #                 color='black', 
    #                 alpha = 0.1, 
    #                 size = 10, 
    #                 legend = False,
    #                ax = ax)    
    
    cutoff_y = np.percentile(ENST_phylo_P['PhyloP'], 1)
    y_range = max(ENST_phylo_P["PhyloP"]) + 1 - cutoff_y
    
    h = ax.hist2d(data=ENST_phylo_P, 
                    x="cds_nt", 
                    y='PhyloP', cmap  = 'binary', bins = [len(smoothed_df) // 25, int(y_range) * 5], alpha =1)
    #fig.colorbar(h[3], ax=ax)
    
    sns.lineplot(data=smoothed_df, x="cds_nt", y='smoothed', color=sns.color_palette('colorblind')[0], ax = ax, lw = 1.5, alpha = 0.7)   


    new_DBD_coords = DBD_coords
    new_AD_coords = AD_coords
    
    # new_DBD_coords = []
    # new_AD_coords = []

    # for start, end in DBD_coords:
    #     new_DBD_coords.append((start * 3, end * 3))
    
    # for start, end in AD_coords:
    #     new_AD_coords.append((start * 3, end * 3))
        
    plot_domain_rectangles(new_DBD_coords, new_AD_coords, RD_coords, max(ENST_phylo_P["PhyloP"]) - min(ENST_phylo_P["PhyloP"]) + 2, 
                           RD = False, AD_text = False, alpha = 0.15,
                          DBD_y = min(ENST_phylo_P["PhyloP"]) - 1, 
                           AD_y = min(ENST_phylo_P["PhyloP"]) - 1, 
                           RD_y = min(ENST_phylo_P["PhyloP"]) - 1, 
                           ax = ax)
    
    sns.despine()
    ax.set_xlim(0, max(smoothed_df["cds_nt"]) + 1)


    ax.set_ylim(cutoff_y, max(ENST_phylo_P["PhyloP"]) + 1)
    
    if alignment == "zoonomia":
        alignment = "Zoonomia"
    elif alignment == "100_verteb":
        alignment = "100 Vertebrates"
    ax.set_ylabel(f"{alignment}\nCDS PhyloP")
    # ax.set_xlabel('CDS Nucloetide', color = sns.color_palette('tab10')[2], labelpad = 0) 

    #ax.set_xticks(np.arange(0, max(ENST_phylo_P["cds_nt"]) + 1, x_axis_spacing))
    #ax.tick_params(axis='x', colors=sns.color_palette('tab10')[2])
    #ax.spines["bottom"].set_color(sns.color_palette('tab10')[2])
    
    ax.axhline(0, color = "gray", linestyle = "--", alpha = 0.7)
    #plt.title(gene)

def plot_prot_trace(gene, uniprotID, ax, x_axis_spacing = 600):
    row = TFs_tbl[TFs_tbl["uniprotID"].str.contains(uniprotID)].iloc[0]
    ENST = row["ENST"]
    DBD_coords = row["DBD_ranges"]
    AD_coords = row["AD_ranges"]
    RD_coords = row["RD_ranges"]

    if "../data/zoonomia_toga_mca/prot_alignment_percent_identities_for_vis/" + gene in glob.glob("../data/zoonomia_toga_mca/prot_alignment_percent_identities_for_vis/*"):
        percent_identities = pd.read_csv("../data/zoonomia_toga_mca/prot_alignment_percent_identities_for_vis/" + gene, index_col = 0)
        #display(percent_identities)
        
        # sns.scatterplot(data = percent_identities, x = "pos", y = "percent_identity", ax = ax, edgecolor = "none", s = 1)
        # display(percent_identities)
    
    
        smoothed = lowess(percent_identities['percent_identity'], 
                          np.arange(len(percent_identities)), 
                          frac=0.01)  # Adjust frac for smoothness
        smoothed_df = pd.DataFrame(smoothed, columns=["pos", 'smoothed'])
    
        # sns.scatterplot(data=percent_identities, 
        #                 x=np.arange(len(percent_identities)), 
        #                 y='percent_identity', 
        #                 color='black', 
        #                 alpha = 0.3, 
        #                 s = 10 * (500 / len(smoothed_df)), 
        #                 legend = False,
        #                ax = ax)  
    
        y_range = 101 - min(min(percent_identities["percent_identity"]),
                        min(smoothed_df['smoothed'])) + 1
    
        
        h = ax.hist2d(data=percent_identities, 
                        x=np.arange(len(percent_identities)), 
                        y='percent_identity', cmap  = 'binary', bins = [len(percent_identities) // 10, 20], alpha = 1)
    
        
        sns.lineplot(data=smoothed_df, x="pos", y='smoothed', color=sns.color_palette('colorblind')[0], ax = ax, lw = 1.5, alpha = 0.7)
    
        sns.despine()
        ax.set_ylabel("Zoonomia\nPercent\nIdentity")
        ax.set_ylim(min(min(percent_identities["percent_identity"]),
                        min(smoothed_df['smoothed'])) - 1, 101)
        ax.set_xlim(0, len(smoothed_df) + 1)
    
        plot_domain_rectangles(DBD_coords, AD_coords, RD_coords, 100, 
                               RD = False, AD_text = False, alpha = 0.15,
                              DBD_y = 0, 
                               AD_y = 0, 
                               RD_y = 0, 
                               ax = ax)
    
        #fig.colorbar(h[3], ax=ax)
        #plt.title(gene)

def plot_mtr_trace(uniprotID, ax, x_axis_spacing = 600):
    row = TFs_tbl[TFs_tbl["uniprotID"].str.contains(uniprotID)].iloc[0]
    #(row)

    ENST = row["ENST"]
    #(ENST)

    DBD_coords = row["DBD_ranges"]
    AD_coords = row["AD_ranges"]
    RD_coords = row["RD_ranges"]
    
    # domain_mtr = []

    # ENSTs = glob.glob(f"../soto_analysis/outputs/mutations/domains_mtr/{ENST}*")

    # for ENST in ENSTs:
    #     try:
    #         df = pd.read_csv(ENST, sep = "\t", header = None)
    #         domain_mtr.append(df)
    #     except pd.errors.EmptyDataError:
    #         print(ENST, " is empty and has been skipped.")

    # domain_mtr = pd.concat(domain_mtr)
    # ad_mtr = domain_mtr[domain_mtr[3] == "AD"]

    cds_mtr = []

    cds_paths = glob.glob(f"../soto_analysis/outputs/mutations/cds_mtr/{ENST}*")
    #print(cds_paths)

    for path in cds_paths:
        file_ENST = path.split("/")[-1].split(".bed")[0]
        if ENST == file_ENST:
            try:
                df = pd.read_csv(path, sep = "\t", header = None)
                cds_mtr.append(df)
            except pd.errors.EmptyDataError:
                print(file_ENST, " is empty and has been skipped.")

    cds_mtr = pd.concat(cds_mtr)
    cds_mtr = cds_mtr.rename(columns = {3: "ENST", 8: "MTR"})    

    ENST_mtr = cds_mtr[cds_mtr["ENST"] == ENST]

    if len(ENST_mtr) > 0:
        if ENST_mtr[4].iloc[0] == -1:
            ascending = False
        else: 
            ascending = True
            
        ENST_mtr = ENST_mtr.sort_values(by = 1, ascending = ascending)
    
        cds_bed = pd.read_csv("../soto_analysis/outputs/mutations/cds_bed_format/" + ENST, sep = "\t", header = None)
        cds_bed["pos"] = [np.arange(s + 1, e + 1) for s, e in zip(cds_bed[1], cds_bed[2])]
    
        full_cds_positions = cds_bed[["pos"]].explode("pos").reset_index(drop = True)
        full_cds_positions = full_cds_positions.sort_values(by = "pos", ascending = ascending)
        #display(full_cds_positions)
    
        merged = pd.merge(full_cds_positions, ENST_mtr, how = "left", left_on = "pos", right_on = 2)
        merged = merged.drop_duplicates(subset = "pos")
        merged["cds_nt"] = np.arange(len(merged)) / 3
        #display(merged)
    
        # ax.axhline(0.841, color = "gray", linestyle = "--", alpha = 0.7)
        
    
        # sns.lineplot(data=merged, x="cds_nt", y='MTR', legend = False, ax = ax, lw = 1, color = sns.color_palette('colorblind')[9])
            
        # sns.despine(ax = ax)
        ax.set_ylabel("Missense\nTolerance\nRatio", rotation = 0, ha = "right", va = "center", labelpad = 10)
    
        
        plot_domain_rectangles(DBD_coords, AD_coords, RD_coords, max(merged["MTR"].dropna()), 
                               RD = False, AD_text = False, alpha = 0.15,
                              DBD_y = 0, 
                               AD_y = 0, 
                               RD_y = 0, 
                               ax = ax)
        
    
        display(merged)

        ax.axhline(0.841, color = "gray", linestyle = "--", alpha = 0.7)
        threshold_plot(data=merged, x="cds_nt", y='MTR', legend = False, ax = ax, lw = 1, 
                       color_below = "silver",
                       color_above = sns.color_palette('colorblind')[9], 
                       threshold = 0.841)
    
    
        min_y = min(merged["MTR"].dropna()) - 0.05
        max_y = max(merged["MTR"].dropna()) + 0.05
        ax.set_ylim(min_y, max_y)
        ax.set_xlim(0, max(merged["cds_nt"]))
    
        center_pos = (0.841 - min_y) / (max_y - min_y)
        add_custom_legend_arrows(top_color = sns.color_palette('colorblind')[-1],
                             top_label = "Non-constrained",
                             bottom_color = sns.color_palette('colorblind')[-3],
                             bottom_label = "Constrained", ax = ax, 
                     center = center_pos, length = 0.45 / 1.5, padding = 0.005)
    
        sns.despine(ax = ax)

import os
def save_structure(uniprotID):
    alphafold_ID = 'AF-' + uniprotID + '-F1'
    database_version = "v6"
    model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_ID}-model_{database_version}.pdb'
    os.system(f'curl {model_url} -o ../data/SFARI_TF_AF/{alphafold_ID}.pdb')


def plot_af_domains(uniprotID, ax, x_axis_spacing = 600):
    save_structure(uniprotID)

    row = TFs_tbl[TFs_tbl["uniprotID"].str.contains(uniprotID)].iloc[0]
    ENST = row["ENST"]

    DBD_coords = row["DBD_ranges"]
    AD_coords = row["AD_ranges"]
    RD_coords = row["RD_ranges"]
    
    p = PDBParser()
    structure = p.get_structure(uniprotID, "../data/SFARI_TF_AF/AF-" + uniprotID + "-F1.pdb")
    try:
        model = structure[0]
        dssp = DSSP(model, "../data/SFARI_TF_AF/AF-" + uniprotID + "-F1.pdb")
    
        pos = []
        code = []
        
        for i in range(len(list(dssp.keys()))):
            pos.append(i)
            a_key = list(dssp.keys())[i]
            code.append(dssp[a_key][2])
    
        codes = pd.DataFrame({"pos" : pos, "code" : code})
        
    
        color_palette = {# Helices
                        "G" : sns.color_palette('colorblind')[0],
                        "H" : sns.color_palette('colorblind')[0],
                        "I" : sns.color_palette('colorblind')[0],
        
                        # Sheets
                        "E" : sns.color_palette('colorblind')[1],
                        "B" : sns.color_palette('colorblind')[1],
        
                        # Poly proline helix
                        "P" : sns.color_palette('colorblind')[2],
        
                        # Other
                        "T" : "none",
                        "S" : "none",
                        "-" : "none"}
    
        #print(DBD_coords)
        
        plot_domain_rectangles(DBD_coords, AD_coords, RD_coords, 1, 
                               RD = False, AD_text = False, alpha = 0.15,
                              DBD_y = 0, 
                               AD_y = 0, 
                               RD_y = 0, 
                               ax = ax)
        
        
        for code in set(codes["code"]):
            intervals =  return_intervals(codes[codes["code"] == code]["pos"])
            for s, e in intervals:
                #print(i)
                color = color_palette[code]
                plot_rectangle(s, e, color, 1, ax, alpha = 1)
                ax.set_xlim(1, len(pos))       
        i+=1
    
            
    
        sns.despine(left = True, ax = ax)
        ax.set_yticks([])
        ax.set_yticklabels([])
    except:
        print("No alphafold structure found")