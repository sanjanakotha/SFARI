import os

d_full_mutations = "../outputs/mutations/cds_clinvar_snv_classified"
d_domain_mutations = "../outputs/mutations/domains_clinvar_snv_classified"

files = os.listdir(d_full_mutations)

with open("../outputs/TFs_table_proteins.txt") as table, open("../outputs/mutations/raw_clinvar_mutations_pathogenic_snv_classified2.txt", "w") as o1:
	table.readline()
	tf_dbd = {}
	tf_ad = {}
	tf_rd = {}
	tf_bif = {}
	tf_full = {}
	tf_effector = {}
	uniprot_tf = {}
	tf_family = {}
	tf_class = {}
	for line in table:
		line = line.rstrip().split("\t")

		tf, dbd, ad, rd, bif, full = line[3], line[6].split(","), line[7].split(","), line[8].split(","), line[9].split(","), line[-1]

		tf_family[line[0]] = line[1]
		uniprot_tf[line[3]] = line[0]
		tf_class[line[0]] = line[2]

		dbd_length, ad_length, rd_length, bif_length, full_length = "NA", "NA", "NA", "NA", "NA"

		full_length = len(line[-1])

		if dbd != ["NA-NA"]:
			dbd_length = 0
			for coord in dbd:
				coord = coord.split("-")
				start, end = int(coord[0]), int(coord[1])
				dbd_length += end - start +1

		if ad != [""]:
			ad_length = 0
			for coord in ad:
				coord = coord.split("-")
				start, end = int(coord[0]), int(coord[1])
				ad_length += end - start + 1
		if rd != [""]:
			rd_length = 0
			for coord in rd:
				coord = coord.split("-")
				start, end = int(coord[0]), int(coord[1])
				rd_length += end - start + 1
		if bif != [""]:
			bif_length = 0
			for coord in bif:
				coord = coord.split("-")
				start, end = int(coord[0]), int(coord[1])
				bif_length += end - start + 1

		effector_length = 0
		if ad_length != "NA":
			effector_length += ad_length
		if rd_length != "NA":
			effector_length += rd_length
		if bif_length != "NA":
			effector_length += bif_length

		tf_full[tf] = full_length * 3

		if ad_length == "NA":
			tf_ad[tf] = ad_length
		else:
			tf_ad[tf] = ad_length * 3

		if rd_length == "NA":
			tf_rd[tf] = rd_length
		else:
			tf_rd[tf] = rd_length * 3

		if bif_length == "NA":
			tf_bif[tf] = bif_length
		else:
			tf_bif[tf] = bif_length * 3

		if dbd_length == "NA":
			tf_dbd[tf] = dbd_length
		else:
			tf_dbd[tf] = dbd_length * 3

		tf_effector[tf] = effector_length *3

	o1.write("TF\tFamily\tClass\tUniprot\tfull\tdbd\tad\trd\tbif\teffector\tfull_length\t
    dbd_length\tad_length\trd_length\tbif_length\teffector_length\tfull_prop\tdbd_prop\t
    ad_prop\trd_prop\tbif_prop\teffector_prop\n")

	for file in files:
		with open(d_full_mutations + "/" + file) as f1, open(d_domain_mutations + "/" + file) as f2:

			name = file.split("_")[0]
			full = []
			for line in f1:
				line = line.rstrip().split("\t")
				if line[-1] != "No-Syn":
					continue
				if "pathogenic" in line[-4].lower():
					full.append((line[2], line[13], line[14]))
			full = len(set(full))
			dbd, ad, rd, bif = [], [], [], []
			for line in f2:
				line = line.rstrip().split("\t")
				if line[-1] != "No-Syn":
					continue
				if line[3] == "DBD" and "pathogenic" in line[-4].lower():
					dbd.append((line[2], line[21], line[22]))

				elif line[3] == "AD" and "pathogenic" in line[-4].lower():
					ad.append((line[2], line[21], line[22]))

				elif line[3] == "RD" and "pathogenic" in line[-4].lower():
					rd.append((line[2], line[21], line[22]))

				elif line[3] == "Bif" and "pathogenic" in line[-4].lower():
					bif.append((line[2], line[21], line[22]))

			dbd = len(set(dbd))
			ad = len(set(ad))
			rd = len(set(rd))
			bif = len(set(bif))

			effector = ad + rd +bif

			full_prop = full/tf_full[name]
			effector_prop = effector/ tf_effector[name]

			if tf_dbd[name] == "NA":
				dbd_prop = "NA"
			else:
				dbd_prop = dbd/tf_dbd[name]

			if tf_ad[name] == "NA":
				ad_prop = "NA"
			else:
				ad_prop = ad/tf_ad[name]

			if tf_rd[name] == "NA":
				rd_prop = "NA"
			else:
				rd_prop = rd/tf_rd[name]

			if tf_bif[name] == "NA":
				bif_prop = "NA"
			else:
				bif_prop = bif/tf_bif[name]

			o1.write("\t".join([uniprot_tf[name], tf_family[uniprot_tf[name]], tf_class[uniprot_tf[name]], name, str(full), str(dbd), str(ad), str(rd), str(bif), str(effector), str(tf_full[name]), str(tf_dbd[name]), str(tf_ad[name]), str(tf_rd[name]), str(tf_bif[name]), str(tf_effector[name]), str(full_prop), str(dbd_prop), str(ad_prop), str(rd_prop), str(bif_prop), str(effector_prop)]) + "\n")
			print(name, dbd, ad, rd, bif)


