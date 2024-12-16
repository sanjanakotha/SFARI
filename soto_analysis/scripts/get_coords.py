import subprocess
with open("../outputs/TFs_table_proteins.txt") as f1:
	f1.readline()
	tfs_n = {}
	data = []
	for line in f1:
		line = line.rstrip().split("\t")
		data.append(line)
		tf, uniprot = line[0], line[3]
		n_dom = 0
		ad_dom, rd_dom, bif_dom = line[7].split(","), line[8].split(","), line[9].split(",")
		if ad_dom != [""]:
			n_dom += len(ad_dom)
		if rd_dom != [""]:
			n_dom += len(rd_dom)
		if bif_dom != [""]:
			n_dom += len(bif_dom)
		if tf not in tfs_n:
			tfs_n[tf] = [uniprot, n_dom]
		else:
			if tfs_n[tf][1] <= n_dom:
				tfs_n[tf] = [uniprot, n_dom]
	i = 1
	for line in data:
		tf, uniprot_id, enst = line[0], line[3], line[5]
		if tf != "ZFP57":
			continue
		print(i, tf, enst)
		command_line = 'grep "' + enst + '" ../raw_files/gencode.v36.annotation.gtf | grep -P "\tCDS\t" > ../outputs/mutations/CDS/' + uniprot_id + '.txt'
		subprocess.call(command_line, shell=True)
		i += 1
