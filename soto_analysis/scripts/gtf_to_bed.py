import subprocess
import os

directory = "../outputs/mutations/CDS"

files = os.listdir(directory)

files2 = os.listdir("../outputs/mutations/cds_bed_format")

for file in files:
	print(file)
	name = file.split(".")[0]
	if name + ".bed" in files2:
		continue
	subprocess.call("gtf2bed < " + directory + "/" + file + " > ../outputs/mutations/cds_bed_format/" + name + ".bed", shell = True)
	subprocess.call('sed -e "s/chr//g" ../outputs/mutations/cds_bed_format/' + name + '.bed > tmpfile && mv tmpfile ../outputs/mutations/cds_bed_format/' + name + '.bed', shell = True)
