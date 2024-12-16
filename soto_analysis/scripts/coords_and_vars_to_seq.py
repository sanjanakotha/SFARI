from Bio.Seq import Seq
from Bio import SeqIO
#from Bio.Alphabet import IUPAC
import sys
# import pyensembl
import os
#https://towardsdatascience.com/a-simple-guide-to-command-line-arguments-with-argparse-6824c30ab1c3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--variants_filename', type=str, required=True)
args = parser.parse_args()

variants_filename = args.variants_filename
variants_filename_short = variants_filename.split(".")[0]

# +
directory = "../outputs/mutations/cds_" + variants_filename_short
files = os.listdir(directory)
o_directory = "../outputs/mutations/cds_" + variants_filename_short + "_snv_classified"

if not os.path.exists(o_directory):
    os.makedirs(o_directory)
        
# -

reverse = {"C":"G", "A":"T", "T":"A", "G":"C"}



# SK: Dictionary of protein names to sequences
proteins = {}
for record in SeqIO.parse("../raw_files/gencode.v36.pc_translations.fa", "fasta"):
    name = record.id.split("|")[1].split(".")[0]
    proteins[name] = str(record.seq)

# +
# SK: Dictionary of protein names to CDS dna transcript

dna_transcripts = {}
for record in SeqIO.parse("../raw_files/gencode.v36.pc_transcripts.fa", "fasta"):
	name = record.id.split("|")[0].split(".")[0]
	record_c = record.id.split("|")
	for i in record_c:
		if "CDS" in i:
			coords = i.replace("CDS:","")
	start = int(coords.split("-")[0])
	end = int(coords.split("-")[1])
	dna_seq = str(record.seq)[start-1:end]
	dna_transcripts[name] = dna_seq
# -

count = 1
for file in files:
    print(count, file)
    count += 1
    with open(directory +  "/" + file) as f1, open(o_directory + "/" + file, "w") as o1:
        
        print_count = 0 # SK
        
        # SK: Reading in the variants from the bed file
        for line in f1:
            # SK: Read in strand and mutated nucleotide
            line = line.rstrip().split("\t")
            strand = line[5]
                
            # mut_nt = line[-4]
            mut_nt = line[-2]
            
            # SK: Check that SNP
            if len(mut_nt) != 1:
                continue
            if len(line[-3]) != 1:
                continue
            if mut_nt == ".":
                continue

            # SK: Not sure if this should be zero or one indexed, assuming one
            # SK 10/25: running into an indexing issue.. not sure why 
            mut_pos = int(line[-4])
            # note = line[9].split(";")[1].split('"')
            # note = line[-1]
            enst = line[4]
            # name = file.split("_")[0]
            name = line[4]
            
            # SK: Reading in the sorted cds genomic coordinates
            list_coords = []
            with open("../outputs/mutations/cds_bed_format/sorted/" + name + ".bed") as f2:
                
                # SK: Iterate through exons in bed file, append tuples of coords to list_coords
                for linex in f2:
                    linex = linex.rstrip().split("\t")
                    #print(linex)
                    start, end = int(linex[1]), int(linex[2])
                    #print(start)
                    #print(end)
                    list_coords.append((start+1, end)) # SK: Bed file -> One-indexed coordinates now
                    
            # SK: Reverse list_coords if on negative strand and reverse translate mutant nucleotide
            if strand == "-":
                list_coords.reverse()
                mut_nt = reverse[mut_nt]
            
            # SK: Get DNA transcript and protein sequence of protein
            nt_seq, translation = dna_transcripts[enst], proteins[enst]

            # Create a list of coords
            total_coords = []
            # SK: Append all individual coordinates in order to total_coords
            for start, end in list_coords:
                if strand == "-":
                    total_coords += list(range(end, start - 1, -1))
         
                else:
                    total_coords += list(range(start, end + 1))

            # SK: dictionary of ith total coordinate to the ith nt_seq
            pos_nt = {} # 1500: A-T-C-G
            for i in range(len(total_coords)):
                if i >= len(nt_seq):
                    print(i)
                    print(nt_seq)
                    print(total_coords[i])
                    print(list_coords)
                pos_nt[total_coords[i]] = nt_seq[i]

            # SK: Breaking total nt coordinates into total codon coordinates
            start = 0
            aa_codon = {} # 300 : [500, 520, 521]
            for i in range(int(len(total_coords)/3)):
                aa_codon[i+1] = total_coords[start:start+3]
                start += 3

            # SK: 
            codon = []
            for k, v in aa_codon.items():
                if mut_pos in v:
                    aa_pos = k
                    codon = v
                    break
                    
            # WT codon
            wt_dict = {} # pos1:nt1, pos2:nt2, pos3:nt3
            for pos in codon:
                wt_dict[pos] = pos_nt[pos]

            # Mutated codon
            mt_dict = dict(wt_dict)
            mt_dict[mut_pos] = mut_nt # pos1:nt1, pos2:nt2, pos3:nt3

            #print(wt_dict, mt_dict)

            # To String
            wt_str = ""
            for nt in wt_dict.values():
                wt_str += nt
            mt_str = ""
            for nt in mt_dict.values():
                mt_str += nt


            wt_aa, mt_aa = Seq(wt_str).translate(), Seq(mt_str).translate()
            line.append(str(wt_aa))
            line.append(str(mt_aa))
            if wt_aa == mt_aa:
                line.append("Syn")
            elif wt_aa != mt_aa and mt_aa != "*":
                line.append("No-Syn")
            elif wt_aa != mt_aa and mt_aa == "*":
                line.append("NoSense")
            o1.write("\t".join(line) + "\n")




