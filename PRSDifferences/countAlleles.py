#!/usr/bin/env python3

def revc(base):
	d = {"A" : "T", "T" : "A", "C" : "G", "G" : "C"}
	return d[base]
import json, os, sys
o = open(sys.argv[4], "w")
o.write("rsid\ttrait\teffect_allele\tmajor_allele\tminor_allele\tCLM_major\tCLM_minor\tCHO_major\tCHO_minor\n")
if len(sys.argv) < 4:
	print(f"countAlleles.py GWAS-catalog-file  clm-freqs-base cho-freqs-base")
s = {"clm" : {}, "cho" : {}}
gwas = {}
with open(sys.argv[1], 'r') as gwas_fh:
	for line in gwas_fh:
		line = line.rstrip().split('\t')
		if len(line) != 3: continue
		elif line[2].upper() not in ["A", "T", "C", "G"]: continue
		elif line[1] not in gwas:
			gwas[line[1]] = {}
			# print(line[1])
			gwas[line[1]]['trait'] = [ line[0] ]
			gwas[line[1]]["allele"] = [ line[2].upper() ]
		else:
			gwas[line[1]]['trait'].append(f"{line[0]}")
			gwas[line[1]]['allele'].append(f"{line[2].upper()}")
print(f"Unique snps in {sys.argv[1]}: {len(gwas.keys())}")
cnt_clm = 0
cnt_cho = 0
for i in range(1,23):
	f1 = sys.argv[2].replace("#",str(i))
	f2 = sys.argv[3].replace("#",str(i))
	with open(f1, 'r') as clm:
		clm.readline()
		for l in clm:
			l = l.rstrip().split()
			if l[1] not in s["clm"] and l[1] in gwas:
				# print(gwas[l[1]])
				cnt_clm += 1
				s["clm"][l[1]] = {}
				s["clm"][l[1]]["a1"] = l[2]
				s["clm"][l[1]]["a2"] = l[3]
				s["clm"][l[1]]["c1"] = int(l[4])
				s["clm"][l[1]]["c2"] = int(l[5])
			else:
				s["clm"][l[1]] = "skip"
	print(f"After filtering chrmosome {i}, CLM has {cnt_clm} SNPs")
	with open(f2, 'r') as cho:
		cho.readline()
		for l in cho:
			l = l.rstrip().split()
			if l[1] not in s["cho"] and l[1] in gwas:
				cnt_cho += 1
				s["cho"][l[1]] = {}
				s["cho"][l[1]]["a1"] = l[2]
				s["cho"][l[1]]["a2"] = l[3]
				s["cho"][l[1]]["c1"] = int(l[4])
				s["cho"][l[1]]["c2"] = int(l[5])
			else:
				s["cho"][l[1]] = "skip"
	print(f"After filtering chrmosome {i}, CHO has {cnt_cho} SNPs")
x = dict()
for k in s["clm"]:
	if k in s["cho"] and s["clm"][k] != "skip" and s["cho"][k] != "skip":
		# print(k)
		if s["cho"][k]["c1"] + s["cho"][k]["c2"] > 150 and s["clm"][k]["c1"] + s["clm"][k]["c2"] > 150:
			for i in range(len(gwas[k]['trait'])):
				if gwas[k]['allele'][i] == s["clm"][k]["a1"] or gwas[k]['allele'][i] == s["clm"][k]["a2"]:
					print(k, "No flip needed", gwas[k]['allele'][i], s["clm"][k]["a1"], s["clm"][k]["a2"])
					pass
				elif revc(gwas[k]['allele'][i]) == s["clm"][k]["a1"] or revc(gwas[k]['allele'][i]) == s["clm"][k]["a2"]:
					gwas[k]['allele'][i] = revc(gwas[k]['allele'][i])
					print(k, "Flip needed",revc(gwas[k]['allele'][i]), "->", gwas[k]['allele'][i], s["clm"][k]["a1"], s["clm"][k]["a2"])
				else:
					# print("Invalid allele")
					continue
				if s["clm"][k]["c1"] > s["clm"][k]["c2"]:
					if s["clm"][k]["a1"] == s["cho"][k]["a1"]:
						# x[k] = {"clm" : {"major" : s["clm"][k]["c1"], "minor" : s["clm"][k]["c1"]},
						# 		"cho" :  {"major" : s["cho"][k]["c1"], "minor" : s["cho"][k]["c1"]} }
						o.write(f'{k}\t{gwas[k]["trait"][i]}\t{gwas[k]["allele"][i]}\t{s["clm"][k]["a1"]}\t{s["clm"][k]["a2"]}\t{s["clm"][k]["c1"]}\t{s["clm"][k]["c2"]}\t{s["cho"][k]["c1"]}\t{s["cho"][k]["c2"]}\n')
					else:
						# x[k] = {"clm" : {"major" : s["clm"][k]["c1"], "minor" : s["clm"][k]["c1"]},
						# 		"cho" :  {"major" : s["cho"][k]["c2"], "minor" : s["cho"][k]["c1"]} }
						o.write(f'{k}\t{gwas[k]["trait"][i]}\t{gwas[k]["allele"][i]}\t{s["clm"][k]["a2"]}\t{s["clm"][k]["a1"]}\t{s["clm"][k]["c1"]}\t{s["clm"][k]["c1"]}\t{s["cho"][k]["c2"]}\t{s["cho"][k]["c1"]}\n')
				else:
					if s["clm"][k]["a1"] == s["cho"][k]["a1"]:
						# x[k] = {"clm" : {"major" : s["clm"][k]["c2"], "minor" : s["clm"][k]["c1"]},
						# 		"cho" :  {"major" : s["cho"][k]["c2"], "minor" : s["cho"][k]["c1"]} }
						o.write(f'{k}\t{gwas[k]["trait"][i]}\t{gwas[k]["allele"][i]}\t{s["clm"][k]["a2"]}\t{s["clm"][k]["a1"]}\t{s["clm"][k]["c2"]}\t{s["clm"][k]["c1"]}\t{s["cho"][k]["c2"]}\t{s["cho"][k]["c1"]}\n')
					else:
						# x[k] = {"clm" : {"major" : s["clm"][k]["c1"], "minor" : s["clm"][k]["c2"]},
						# 		"cho" :  {"major" : s["cho"][k]["c1"], "minor" : s["cho"][k]["c2"]} }
						o.write(f'{k}\t{gwas[k]["trait"][i]}\t{gwas[k]["allele"][i]}\t{s["clm"][k]["a1"]}\t{s["clm"][k]["a2"]}\t{s["clm"][k]["c1"]}\t{s["clm"][k]["c2"]}\t{s["cho"][k]["c1"]}\t{s["cho"][k]["c2"]}\n')



# with open(f"{sys.argv[3]}.json", 'w') as j:
# 	json.dump(x,j)
