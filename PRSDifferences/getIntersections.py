#!/bin/env python3
import sys
import os
if len(sys.argv) != 4:
	print(f"{sys.argv[0]} traitFile traitFolder snpFile-Folder outFolder")

if not os.path.exists(sys.argv[3]):
	os.makedirs(sys.argv[3], exist_ok=True)

clm_list, cho_list = [], []

clmFile = f"{sys.argv[3]}/{sys.argv[1]}".replace(".tsv",".clm")
choFile = f"{sys.argv[3]}/{sys.argv[1]}".replace(".tsv",".cho")

with open(clmFile, "r") as CLM:
	CLM.readline()
	for line in CLM:
		rs = line.split("\t")[0]
		clm_list.append(rs)

with open(choFile, "r") as CHO:
	CHO.readline()
	for line in CHO:
		rs = line.split("\t")[0]
		cho_list.append(rs)

print(cho_list)
print(clm_list)

intersection = list(set(cho_list) & set(clm_list))

with open(f"{sys.argv[2]}/{sys.argv[1]}", "r") as traits:
	with open(f"{sys.argv[4]}/{sys.argv[1]}", "w") as updatedTraits:
		for line in traits:
			rs = line.split("\t")[0]
			if rs in intersection:
				updatedTraits.write(line)
