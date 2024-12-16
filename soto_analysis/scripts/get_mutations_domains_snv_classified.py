import subprocess
import os

#https://towardsdatascience.com/a-simple-guide-to-command-line-arguments-with-argparse-6824c30ab1c3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--domain', type=str, required=True)
parser.add_argument('--variants_filename', type=str, required=True)
args = parser.parse_args()

variants_filename = args.variants_filename
variants_filename_short = variants_filename.split(".")[0]


# +
d_mutations = "../outputs/mutations/cds_" + variants_filename_short + "_snv_classified"
folder = args.domain

if folder == "none":
    d_domains = "../outputs/mutations/domains_bed_format"
else:
    d_domains = "../outputs/mutations/domains_bed_format" + "/" + folder
files = os.listdir(d_domains)
# -


#SK 
import pandas as pd
uniprotID_ENST_mapping = pd.read_csv("../../output/TFs_with_ENST.csv")#data/SFARI_TFs_with_ENST_corrected.csv")
uniprotID_ENST_mapping = uniprotID_ENST_mapping[["uniprotID", "ENST"]]
uniprotID_ENST_mapping["ENST"] = uniprotID_ENST_mapping["ENST"].str.split(".").str[0]
uniprotID_ENST_mapping_dict= dict(zip(uniprotID_ENST_mapping["uniprotID"], uniprotID_ENST_mapping["ENST"]))
#uniprotID_ENST_mapping_dict['O60479'] = 'ENST00000434704'
uniprotID_ENST_mapping_dict

i = 1
for file in files:
    name = file.split(".")[0] # gene
    if (name in uniprotID_ENST_mapping_dict.keys()):
        print(i, file) # gene.bed
        ENST = uniprotID_ENST_mapping_dict[name]
        if folder == "none":
            subprocess.call("bedtools intersect -a " + d_domains + "/" + file + " -b " + d_mutations + "/" + ENST + ".bed" + " -wb > ../outputs/mutations/domains_" + variants_filename_short + "_snv_classified/" + ENST + ".bed", shell = True)
        else:
            subprocess.call("bedtools intersect -a " + d_domains + "/" + file + " -b " + d_mutations + "/" + ENST + ".bed" + " -wb > ../outputs/mutations/domains_" + variants_filename_short + "_snv_classified/" + folder + "/" + ENST + ".bed", shell = True)
        i+=1
    else:
        print(name + " does not have an ENST code")


