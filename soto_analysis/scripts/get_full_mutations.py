import subprocess
import os

directory = "../outputs/mutations/cds_bed_format"

#https://towardsdatascience.com/a-simple-guide-to-command-line-arguments-with-argparse-6824c30ab1c3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--variants_filename', type=str, required=True)
parser.add_argument('--sort_TFs', type=str)
parser.add_argument('--sort_variants', type=str)
parser.add_argument('--variants_folder', type = str)
args = parser.parse_args()

variants_filename = args.variants_filename
if args.variants_folder:
    variants_folder = args.variants_folder
else:
    variants_folder = "../raw_files/"

variants_filepath = variants_folder + variants_filename 
variants_filename_short = variants_filename.split(".")[0]
files = os.listdir(directory)

### SK: Sorting the TF bed files
if args.sort_TFs:
    print("Sorting TF bed files.")
    i = 1
    for file in files:
        print(i, file) # gene.bed
        name = file.split(".")[0] # gene
        if not os.path.exists("../outputs/mutations/cds_bed_format/sorted/"):
            os.makedirs("../outputs/mutations/cds_bed_format/sorted/")
        # print("bedtools sort -i " + directory + "/" + file + " > ../outputs/mutations/cds_bed_format_sorted/" + name + ".bed")
        subprocess.call("bedtools sort -i " + directory + "/" + file + " > ../outputs/mutations/cds_bed_format/sorted/" + name + ".bed", shell = True)
        i+=1
    
    print("Done sorting the TF bed files!")

# +
# Sorting the variants
if args.sort_variants:
    print("Sorting variants. at " + variants_filepath)
    if not os.path.exists(variants_folder + "sorted/"):
        os.makedirs(variants_folder + "sorted/")
    subprocess.call("bedtools sort -i " + variants_filepath + " > " + variants_folder + "sorted/" + variants_filename , shell = True)
    print("Done sorting variants!")

variants_filepath = variants_folder + "sorted/" + variants_filename 



# +
print(variants_filepath)

i = 1
print("Starting intersecting CDS bed files and variants.")
if not os.path.exists("../outputs/mutations/cds_" + variants_filename_short):
    os.makedirs("../outputs/mutations/cds_" + variants_filename_short)

for file in files:
    path = os.path.join(directory, file)
    if os.path.isdir(path):
        # skip directories
        print("skipping directory! " + str(path))
        continue
        

    print(i, file) # gene.bed
    
    #print(directory + "/sorted/" + file)
    #print(variants_filepath)
    name = file.split(".")[0] # gene

    
    subprocess.call("bedtools intersect -a " + directory + "/sorted/" + file + ".bed -b " + variants_filepath + " -wb -sorted -nonamecheck > ../outputs/mutations/cds_" + variants_filename_short + "/" + name + ".bed", shell = True)
    i+=1
print("Done intersecting.")
# -


