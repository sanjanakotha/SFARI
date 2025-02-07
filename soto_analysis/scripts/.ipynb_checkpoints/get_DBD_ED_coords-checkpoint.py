import math
from itertools import groupby
from operator import itemgetter
import os

def split_coords(data, strand):
	# https://stackoverflow.com/questions/3149440/splitting-list-based-on-missing-numbers-in-a-sequence
	output = []
	if strand == "-":
		data.reverse()
	for k, g in groupby(enumerate(data), lambda i_x: i_x[0] - i_x[1]):
		output.append(list(map(itemgetter(1), g)))
	return output

directory = "../outputs/mutations/cds_bed_format"
files = os.listdir(directory)

import sys
with open("../outputs/all_TFs_table_proteins.txt") as f1:
    f1.readline()
    for line in f1:
        line = line.rstrip().split("\t")
        # print(line)
        tf, dbd_coords, ad_coords, rd_coords, bif_coords, length = line[3], line[6].split(","), line[7].split(","), line[8].split(","), line[9].split(","), len(line[-1])
        
        if tf == "Q15583-2":
            continue
        ensg, enst = line[4], line[5]
        
        # Get coordinates of domains
        dbd_dom, ad_dom, rd_dom, bif_dom = [], [], [], []

        if dbd_coords != ["NA-NA"]:
            for coord in dbd_coords:
                coord = coord.split("-")
                start, end = int(coord[0]), int(coord[1])
                dbd_dom.append((start, end))
        if ad_coords != [""]:
            for coord in ad_coords:
                coord = coord.split("-")
                start, end = int(coord[0]), int(coord[1])
                ad_dom.append((start, end))
        if rd_coords != [""]:
            for coord in rd_coords:
                coord = coord.split("-")
                start, end = int(coord[0]), int(coord[1])
                rd_dom.append((start, end))
        if bif_coords != [""]:
            for coord in bif_coords:
                coord = coord.split("-")
                start, end = int(coord[0]), int(coord[1])
                bif_dom.append((start, end))

        all_doms = {"DBD": dbd_dom, "AD": ad_dom, "RD": rd_dom, "Bif": bif_dom}

        
        with open("../outputs/mutations/cds_bed_format/" + enst) as f2:

            coords = []
            nt_length = 0
            for line in f2:
                line = line.rstrip().split("\t")
                print(line)
                strand = line[5]
                chrom = line[0]
                # If the negative strand: switch the end and start of the exons  
                if strand == "-":
                    start, end = int(line[2]), int(line[1]) + 1
                else:
                    start, end = int(line[1]) + 1, int(line[2])

                coords.append((start, end))
                nt_length += abs(end - start) + 1

            aa_length = int(nt_length/3)
            #if length != aa_length:
                #print(tf, length, aa_length)
                #continue

            n_exons = len(coords)
            #SK: commented next two lines out
#             if strand == "-":
#                 coords.reverse()
            # Sort in increasing order 
            if strand == "-":
                coords.sort(key=lambda _: -1 * _[0])


            #print(tf, n_exons, strand)

            dict_aa_coords = {}
            all_pos = []
            for start, end in coords:
                if strand == "+":
                    all_pos += list(range(start, end + 1))
                elif strand == "-":
                    all_pos += list(range(start, end -1, -1))
            start_codon = 0
            for i in range(aa_length):
                end_codon = start_codon +3
                dict_aa_coords[i+1] = all_pos[start_codon : end_codon]
                start_codon += 3
                
            with open("../outputs/mutations/domains_bed_format/" + tf, "w") as o1:
                print(tf)
                for name, domain in all_doms.items():
                    if domain == []:
                        continue
                    #print(name, domain)
                    for start_dom, end_dom in domain:
                        #print(start_dom, end_dom, domain)
                        dom_pos = []
                        for aa in range(start_dom, end_dom +1):
                            dom_pos += dict_aa_coords[aa]
                        for cds_dom in split_coords(dom_pos, strand):
                            start_bed = cds_dom[0] - 1
                            end_bed = cds_dom[-1]
                            o1.write(str(chrom) + "\t" + str(start_bed) + "\t" + str(end_bed) + "\t" + name + "\t" + ensg + "\t.\t" + strand + "\t" + enst + "\n")
                            #print(tf, name, start_bed, end_bed)


