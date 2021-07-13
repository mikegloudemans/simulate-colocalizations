#!/usr/bin/python
# Author: Mike Gloudemans
# Date created: 7/16/2018

import os
import sys
import config
import datetime
from shutil import copyfile
from io import StringIO
import gzip
import operator
import random
import math
import subprocess
import pandas as pd
import numpy as np
from scipy import stats

config_file = sys.argv[1]

locus_num = sys.argv[2]

settings = config.load_config(config_file)

base_output_dir = f"output/simulations/{settings['out_dir_group']}/"
os.makedirs(f"{base_output_dir}/hg19/gwas", exist_ok=True)
os.makedirs(f"{base_output_dir}/hg19/eqtl", exist_ok=True)
os.makedirs(f"{base_output_dir}/hg38/gwas", exist_ok=True)
os.makedirs(f"{base_output_dir}/hg38/eqtl", exist_ok=True)

copyfile(config_file, f"{base_output_dir}/settings_used.config")

tmp_dir = f"tmp/{locus_num}"
os.makedirs(tmp_dir, exist_ok=True)

def main():

	# NOTE: There is a possible race condition here; try to fix that to 
	# make sure we're not appending to the file before the header's written
	if locus_num == "1":
		with open(f"{base_output_dir}/answer_key.txt", "w") as w:
			w.write("test_case\tseed_chrom\tseed_pos\tcausal_variants\tcases_n\tcontrols_n\teqtl_n\tsnps_n\n")

	# Get possible GWAS loci upfront, so we don't
	# waste any more time on this

	locus_bank = get_possible_loci(settings)

	random.seed(int(locus_num))
	
	# NOTE: This could be dangerous because it could lead to infinite loops

	# Loop until we get a variant that works
	while True:

		# Store settings for current run, for easy reference
		settings["current_run"] = {}

		print("Starting new locus...")
		# Choose a locus 
		locus = random.choice(locus_bank)
		print(locus)

		# Figure out if there's going to be a causal GWAS variant
		gwas_log_odds = get_gwas_log_odds(settings)

		# Get all haplotypes using HAPGEN2
		gwas_effect_sizes = run_hapgen2(settings, locus, gwas_log_odds)

		print("Got GWAS effect sizes")

		settings["current_run"]["gwas_effect_sizes"] = gwas_effect_sizes
		if gwas_effect_sizes == "Bad variant" or gwas_effect_sizes == "Fail":
			# Make sure HAPGEN2 actually succeeded, i.e. it was a valid site
			print(gwas_effect_sizes)
			continue

		# Get effect sizes, using LD pairings at restrictions
		eqtl_effect_sizes = get_eqtl_effect_sizes(settings, gwas_effect_sizes)

		print("Got eQTL effect sizes")

		settings["current_run"]["eqtl_effect_sizes"] = eqtl_effect_sizes
		if eqtl_effect_sizes == "Impossible LD matrix":
			continue

		eqtl_phenotypes = get_expression_phenotypes(eqtl_effect_sizes)
		print("Got eQTL phenotypes")

		if "simulate_peer" in settings and settings["simulate_peer"] == "True":
			peer_factors = simulate_peer_factors(eqtl_phenotypes, locus_num)

		# Finally, create summary statistics
		sum_stats = make_sum_stats(eqtl_phenotypes)
		print("made sumstats")
		(eqtls, gwas) = sum_stats

		if "simulate_peer" in settings and settings["simulate_peer"] == "True":
			peer_eqtls = call_eqtls_with_peer(eqtl_phenotypes, peer_factors)
			write_peer_sumstats(peer_eqtls, locus_num)

		write_sumstats(eqtls, gwas, locus_num)
		write_answer_key(gwas_effect_sizes, eqtl_effect_sizes, locus, locus_num)
		write_expression_phenotypes(eqtl_phenotypes, locus_num, locus)

		print("wrote expression phenotypes")

		if not "no_genotypes" in settings or settings["no_genotypes"] == "False":
			write_genotypes_as_vcf(locus_num, locus)

		print("wrote genotypes as vcf")

		break

	# liftOver to hg38
	if not "liftover" in settings or settings["liftover"] == "True":
		run_liftover(settings)

def get_eqtl_effect_sizes(settings, gwas_effect_sizes):

	eqtl_effect_sizes = [0] * len(gwas_effect_sizes)
	potential_eqtl_effect = random.uniform(settings["eqtl_min_effect_size"], settings["eqtl_max_effect_size"])
	if max([abs(ges) for ges in gwas_effect_sizes]) > 0:
		# GWAS is causal
		gwas_causal_indices = [i for i,ges in enumerate(gwas_effect_sizes) if abs(ges) > 0]

		if random.random() < settings["p_eqtl_causal_given_gwas_causal"]:
			# eQTL is causal

			if random.random() < settings["p_same_causal_given_both_causal"]:
				# Same causal variant for both
				eqtl_effect_sizes[gwas_causal_indices[0]] = potential_eqtl_effect
			else:
				# Different eQTL causal variant
				causal_eqtl_index = random.randint(len(eqtl_effect_sizes) * 1 // 4, len(eqtl_effect_sizes) * 3 // 4)
				eqtl_effect_sizes[causal_eqtl_index] = potential_eqtl_effect
	else:
		# GWAS is not causal
		if random.random() < settings["p_eqtl_causal_given_no_gwas_causal"]:
			# eQTL is causal
			# Choose causal eQTL; don't let it be right at the edge of the measured region though
			causal_eqtl_index = random.randint(len(gwas_effect_sizes) * 1 // 4, len(gwas_effect_sizes) * 3 // 4)
			eqtl_effect_sizes[causal_eqtl_index] = potential_eqtl_effect

	return eqtl_effect_sizes

# Get the set of loci that are hits in the given GWAS files
def get_possible_loci(settings):
	possible_loci = []
	for gwas_file in settings["gwas_reference_files"]:

		with gzip.open(gwas_file) as f:
			header = f.readline().strip().decode("utf-8").split()

			pval_index = header.index("pvalue")
			chr_index = header.index("chr")
			snp_pos_index = header.index("snp_pos")

			all_snps = []

			for line in f:
				data = line.strip().decode("utf-8").split("\t")
				try:
					pvalue = float(data[pval_index])
				except:
					continue
				chr = data[chr_index]
				pos = int(data[snp_pos_index])
				if pvalue > settings["gwas_threshold"]:
					continue

				all_snps.append((chr, pos, pvalue))
	   
		# For now, include only autosomal SNPs.
		filtered = []
		for s in all_snps:
			if "chr" in str(s[0]):
				filtered.append((s[0][3:], s[1], s[2]))
			else:
				filtered.append((s[0], s[1], s[2]))

		all_snps = sorted(filtered, key=operator.itemgetter(2)) 
		
		# Go through the list of SNPs in order, adding the ones
		# passing our criteria.
		for snp in all_snps:

			# For now, ignore a SNP if it's in the MHC region -- this
			# would require alternative methods.
			if (snp[0] == "6") and snp[1] > 25000000 and snp[1] < 35000000:
					continue

			# Before adding a SNP, make sure it's not right next
			# to another SNP that we've already selected.
			skip = False
			for kept_snp in possible_loci:
					if kept_snp[0] == snp[0] and abs(kept_snp[1] - snp[1]) < settings["window"]:
							skip = True
							break
			if not skip:
				possible_loci.append(snp)

	return possible_loci

def run_hapgen2(settings, locus, gwas_log_odds):

	# Get the data from 1000 Genomes VCF near locus
	filename = f"data/1KG/hg19/ALL.chr{locus[0]}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

	with gzip.open(f"{filename}") as f:
		for line in f:
			clean = line.strip().decode('utf-8')
			if clean.startswith("#CHROM"):
				header = clean.strip().split()
				break
	stream = StringIO(subprocess.run(f"tabix {filename} {locus[0]}:{locus[1] - settings['window']}-{locus[1] + settings['window']}".split(), capture_output=True).stdout.decode("utf-8"))
	vcf = pd.read_csv(stream, sep="\t", names=header)

	# Remove variants with position appearing multiple times
	dup_counts = {}
	for pos in vcf["POS"]:
			dup_counts[pos] = dup_counts.get(pos, 0) + 1
	vcf["dup_counts"] = [dup_counts[pos] for pos in vcf['POS']]
	vcf = vcf[vcf["dup_counts"] == 1]

	# Filter down to a min MAF of 0.05 for now
	def maf(x):
		af = x.split("AF=")[1].split(";")[0]
		if "," in af:
			return 0
		af = float(af)
		maf = min(abs(1-af), abs(af))
		return maf
	vcf['MAF'] = vcf['INFO'].apply(maf)
	vcf = vcf[vcf['MAF'] > 0.05]

	# Make sure the variant of interest hasn't been filtered out;
	# if it has, then get a new one.
	if int(locus[1]) not in set(list(vcf['POS'])):
		return "Bad variant"
	
	# Format it for HAPGEN
	# Write .haps file
	vcf_genos = vcf.iloc[:,10:(vcf.shape[1]-2)]
	with open(f"{tmp_dir}/hapgen2.haps", "w") as w:
		for i, row in vcf_genos.iterrows():
			w.write("|".join(list(row)).replace("|", " ") + "\n")

	# Write .leg file
	with open(f"{tmp_dir}/hapgen2.leg", "w") as w:
		w.write("rs position X0 X1\n")
		for i, row in vcf.iterrows():
			w.write(" ".join([str(r) for r in [row['ID'], row['POS'], row['REF'], row['ALT']]]) + "\n")

	# Write .map centimorgans recombination rate file
	map_file = "data/genetic-map/extended_genetic_map_GRCh37_chr{0}.txt.sorted.gz".format(locus[0])

	#with gzip.open(f"{map_file}") as f:
	#	head = f.readline().decode('utf-8').strip().split()
	
	vcf_pos_set = set(list(vcf['POS']))

	f = StringIO(subprocess.run(f"tabix {map_file} {locus[0]}:{locus[1] - settings['window']}-{locus[1] + settings['window']}".split(), capture_output=True).stdout.decode("utf-8"))
	with open(f"{tmp_dir}/hapgen2.map", "w") as w:
		#with open(genetic_map) as f:
			w.write(f.readline())
			for line in f:
				if line == "":
					continue
				data = line.strip().split()
				if int(data[1]) in vcf_pos_set:
					w.write(line)

	gwas_control_sample_size = int(math.floor(10**random.uniform(settings["gwas_min_control_log_sample_size"], settings["gwas_max_control_log_sample_size"])))
	gwas_case_sample_size = int(math.floor(10**random.uniform(settings["gwas_min_case_log_sample_size"], settings["gwas_max_case_log_sample_size"])))
	eqtl_sample_size = int(math.floor(10**random.uniform(settings["eqtl_min_log_sample_size"], settings["eqtl_max_log_sample_size"])))
	settings["current_run"]["gwas_control_sample_size"] = gwas_control_sample_size
	settings["current_run"]["gwas_case_sample_size"] = gwas_case_sample_size
	settings["current_run"]["eqtl_sample_size"] = eqtl_sample_size

	# Then run HAPGEN2 to get genotypes
	try:
		subprocess.run(f'''bin/hapgen2 -m {tmp_dir}/hapgen2.map -l {tmp_dir}/hapgen2.leg -h {tmp_dir}/hapgen2.haps -o {tmp_dir}/hapgen2_gwas.out -dl {locus[1]} 1 {10 ** gwas_log_odds} {10 ** (gwas_log_odds * 2)} -n {gwas_control_sample_size} {gwas_case_sample_size}'''.split())
		subprocess.run(f'''bin/hapgen2 -m {tmp_dir}/hapgen2.map -l {tmp_dir}/hapgen2.leg -h {tmp_dir}/hapgen2.haps -o {tmp_dir}/hapgen2_eqtl.out -dl {locus[1]} 1 1 1 -n {eqtl_sample_size} 1'''.split())
	except:
		return "Fail"

	# Write VCF to tmp file in case we need to compute LD
	vcf.to_csv(f'{tmp_dir}/tmp.vcf', sep="\t", index=False, header=True)
	settings["current_run"]["rsids"] = vcf['ID']

	gwas_effect_sizes = [0] * vcf.shape[0]
	gwas_effect_sizes[list(vcf['POS']).index(locus[1])] = gwas_log_odds

	return gwas_effect_sizes
	
def get_gwas_log_odds(settings):
  
	if random.random() < settings["p_gwas_causal"]:
		return random.uniform(settings["gwas_min_log_odds"], settings["gwas_max_log_odds"]) * random.choice([-1,1])
	else:
		return 0

def get_expression_phenotypes(eqtl_effect_sizes):

	# Read haplotypes for each individual
	haps = pd.read_csv(f"{tmp_dir}/hapgen2_eqtl.out.controls.haps", sep=" ", header=None)
	genos = np.add(haps.iloc[:, range(0,haps.shape[1]-1,2)].values, haps.iloc[:, range(1,haps.shape[1]-1,2)].values)

	means = np.matmul(genos.T, eqtl_effect_sizes)
	phenos = np.random.normal(means, 1, len(means))

	return phenos

def simulate_peer_factors(phenos, index):

	all_peer_factors = []

	for frac in settings["peer_fraction_var_explained"]:
		pheno_var = np.var(phenos)
		total_var = pheno_var / frac
		peer_factors = np.random.normal(phenos, total_var-pheno_var, len(phenos))

		with open(f"{base_output_dir}/hg19/eqtl/peer_factors{locus_num}_{frac}.txt", "w") as w:
			header = "\t".join(["ID" + str(i) for i in range(len(phenos))]) + "\n"
			w.write(header)
			w.write("\t".join([str(pf) for pf in peer_factors]))
			w.write("\n")

		all_peer_factors.append(peer_factors)

	return all_peer_factors

def write_expression_phenotypes(phenos, index, locus):
	# We do this because RTC needs to re-call eQTLs with various SNPs regressed out
	with open(f"{base_output_dir}/hg19/eqtl/eqtl_phenotypes{locus_num}.bed", "w") as w:

		header = "#chr\tstart\tend\tgene\tlength\tstrand\t" + "\t".join(["ID" + str(i) for i in range(len(phenos))]) + "\n"
		w.write(header)

		# Now put in the phenotypes for our fake gene
		w.write(f"chr{locus[0]}\t{locus[1]}\t{int(locus[1]) + 500}\tnameless_gene\t500\t+\t")
		w.write("\t".join([str(p) for p in phenos]))
		w.write("\n")

	# bgzip and tabix it
	subprocess.run(f"bgzip -f {base_output_dir}/hg19/eqtl/eqtl_phenotypes{locus_num}.bed".split())
	subprocess.run(f"tabix -f -S 1 -s 1 -b 2 -e 3 {base_output_dir}/hg19/eqtl/eqtl_phenotypes{locus_num}.bed.gz".split())

def write_genotypes_as_vcf(index, locus):
	
	metadata = pd.read_csv(f"{tmp_dir}/hapgen2_eqtl.out.controls.gen", sep=" ", header=None)
	haps = pd.read_csv(f"{tmp_dir}/hapgen2_eqtl.out.controls.haps", sep=" ", header=None)
	
	# We do this because RTC needs to re-call eQTLs with various SNPs regressed out
	with open(f"{base_output_dir}/hg19/eqtl/eqtl_genotypes{locus_num}.vcf", "w") as w:

		w.write("##fileformat=VCFv4.1\n") # This is essential; otherwise QTLtools will not work
		header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(["ID" + str(i) for i in range(haps.shape[1] // 2)]) + "\n"
		w.write(header)

		# Now put in our simulated genotypes

		# NOTE: For large sample sizes this can be incredibly slow.
		# Actually though it's way faster now that it's fixed
		'''
		for i in range(metadata.shape[0]):
			w.write(f"chr{locus[0]}\t{metadata.iloc[i,2]}\t{locus[0]}_{metadata.iloc[i,2]}\t{metadata.iloc[i,3]}\t{metadata.iloc[i,4]}\t100\tPASS\tNONE\tGT")
			for j in range(haps.shape[1] // 2):
				w.write("\t" + "|".join([str(int(s)) for s in list(haps.iloc[i,2*j:(2*(j+1))])]))
			w.write("\n")
		'''

		for i in range(metadata.shape[0]):
			metaslice = list(metadata.iloc[i,:])
			w.write(f'''chr{locus[0]}\t{metaslice[2]}\t{locus[0]}_{metaslice[2]}\t{metaslice[3]}\t{metaslice[4]}\t100\tPASS\tNONE\tGT''')
			hapslice = haps.iloc[i,:]
			chunk = ""
			for j in range(haps.shape[1] // 2):
				chunk += "\t" + "|".join([str(int(s)) for s in list(hapslice[2*j:(2*(j+1))])])
			w.write(chunk)
			w.write("\n")


	# bgzip and tabix it
	subprocess.run(f"bgzip -f {base_output_dir}/hg19/eqtl/eqtl_genotypes{locus_num}.vcf".split())
	subprocess.run(f"tabix -f -S 1 -s 1 -b 2 -e 2 {base_output_dir}/hg19/eqtl/eqtl_genotypes{locus_num}.vcf.gz".split())

def make_sum_stats(eqtl_phenotypes):
	eqtls = call_eqtls(eqtl_phenotypes)
	gwas = run_gwas()
	return (eqtls, gwas)

def call_eqtls(eqtl_phenotypes):
	# Load genotypes
	# Read haplotypes for each individual
	haps = pd.read_csv(f"{tmp_dir}/hapgen2_eqtl.out.controls.haps", sep=" ", header=None)
	genos = np.add(haps.iloc[:, range(0,haps.shape[1]-1,2)].values, haps.iloc[:, range(1,haps.shape[1]-1,2)].values)

	# Call eQTLs
	eqtls = []
	for i in range(genos.shape[0]):
		predictors = genos[i,:]
		slope, intercept, r_value, p_value, std_err = stats.linregress(predictors, eqtl_phenotypes)
		if np.isnan(slope):
			slope = 0
			std_err = 1
		eqtls.append((slope, std_err, p_value))
	return eqtls

def call_eqtls_with_peer(eqtl_phenotypes, peer):

	# Load genotypes
	# Read haplotypes for each individual
	haps = pd.read_csv(f"{tmp_dir}/hapgen2_eqtl.out.controls.haps", sep=" ", header=None)
	genos = np.add(haps.iloc[:, range(0,haps.shape[1]-1,2)].values, haps.iloc[:, range(1,haps.shape[1]-1,2)].values)

	all_eqtls = []

	for p in peer:
		# Get residuals for gene expression
		slope, intercept, r_value, p_value, std_err = stats.linregress(p, eqtl_phenotypes)
		preds = p * slope + intercept
		residuals = eqtl_phenotypes - preds

		# Call eQTLs
		eqtls = []
		for i in range(genos.shape[0]):
			predictors = genos[i,:]
			slope, intercept, r_value, p_value, std_err = stats.linregress(predictors, residuals)
			if np.isnan(slope):
				slope = 0
				std_err = 1
			eqtls.append((slope, std_err, p_value))

		all_eqtls.append(eqtls)

	return all_eqtls



def run_gwas():
	# Load genotypes
	# Read haplotypes for each individual
	case_haps = pd.read_csv(f"{tmp_dir}/hapgen2_gwas.out.cases.haps", sep=" ", header=None)
	case_genos = np.add(case_haps.iloc[:, range(0,case_haps.shape[1]-1,2)].values, case_haps.iloc[:, range(1,case_haps.shape[1]-1,2)].values)
	control_haps = pd.read_csv(f"{tmp_dir}/hapgen2_gwas.out.controls.haps", sep=" ", header=None)
	control_genos = np.add(control_haps.iloc[:, range(0,control_haps.shape[1]-1,2)].values, control_haps.iloc[:, range(1,control_haps.shape[1]-1,2)].values)

	# Run GWAS (for now, maybe just do it as a linear regression too? Makes it easily adaptable for continous phenotypes)
	gwas = []
	for i in range(case_genos.shape[0]):
		predictors = list(case_genos[i,:]) + list(control_genos[i,:])
		gwas_response = [1] * case_genos.shape[1] + [0] * control_genos.shape[1]
		slope, intercept, r_value, p_value, std_err = stats.linregress(predictors, gwas_response)
		# Catch issues in cases where one allele never occurred
		if np.isnan(slope):
			slope = 0
			std_err = 1
		gwas.append((slope, std_err, p_value))
	return gwas

def write_sumstats(eqtls, gwas, index):
	#header = subprocess.check_output("cat /users/mgloud/projects/coloc_comparisons/tmp/tmp.vcf 2> /dev/null | grep \\#CHROM", shell=True).strip().split()
	vcf = pd.read_csv(f"{tmp_dir}/tmp.vcf", sep="\t")

	# Filter down to a min MAF of 0.05 for now
	def alt_af(x):
		af = x.split("AF=")[1].split(";")[0]
		if "," in af:
			return 0
		af = float(af)
		return(af)
	vcf['AF'] = vcf['INFO'].apply(alt_af)
 
	with open(f"{base_output_dir}/hg19/gwas/gwas_sumstats{locus_num}.txt", "w") as w:
		w.write("rsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref\talt\teffect_af\tbeta\tse\tzscore\tpvalue\tn_cases\tn_controls\n")
		for i in range(vcf.shape[0]):
			var = vcf.iloc[i, :]
			g = gwas[i]
			variant_id = f'''{var['#CHROM']}_{var["POS"]}_{var["REF"]}_{var["ALT"]}_hg19'''
			w.write(f'''{var['ID']}\t{variant_id}\thg19\t{var['#CHROM']}\t{var["POS"]}\t{var["REF"]}\t{var["ALT"]}\t{var["AF"]}\t{g[0]}\t{g[1]}\t{g[0] / g[1]}\t{g[2]}\t{settings["current_run"]["gwas_case_sample_size"]}\t{settings["current_run"]["gwas_control_sample_size"]}\n''')

	with open(f"{base_output_dir}/hg19/eqtl/eqtl_sumstats{locus_num}.txt", "w") as w:
		w.write("feature\trsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref\talt\teffect_af\tbeta\tse\tzscore\tpvalue\tN\n")
		for i in range(vcf.shape[0]):
			var = vcf.iloc[i, :]
			e = eqtls[i]
			variant_id = f'''{var['#CHROM']}_{var["POS"]}_{var["REF"]}_{var["ALT"]}_hg19'''
			w.write(f'''nameless_gene\t{var['ID']}\t{variant_id}\thg19\t{var['#CHROM']}\t{var["POS"]}\t{var["REF"]}\t{var["ALT"]}\t{var["AF"]}\t{e[0]}\t{e[1]}\t{e[0] / e[1]}\t{e[2]}\t{settings["current_run"]["eqtl_sample_size"]}\n''')

	# bgzip and tabix
	subprocess.run(f"bgzip -f {base_output_dir}/hg19/gwas/gwas_sumstats{locus_num}.txt".split())
	subprocess.run(f"bgzip -f {base_output_dir}/hg19/eqtl/eqtl_sumstats{locus_num}.txt".split())
	subprocess.run(f"tabix -f -S 1 -s 4 -b 5 -e 5 {base_output_dir}/hg19/gwas/gwas_sumstats{locus_num}.txt.gz".split())
	subprocess.run(f"tabix -f -S 1 -s 5 -b 6 -e 6 {base_output_dir}/hg19/eqtl/eqtl_sumstats{locus_num}.txt.gz".split())

def write_peer_sumstats(eqtls, index):
	vcf = pd.read_csv("/users/mgloud/projects/coloc_comparisons/{tmp_dir}/tmp.vcf", sep="\t")

	# Filter down to a min MAF of 0.05 for now
	def alt_af(x):
		af = x.split("AF=")[1].split(";")[0]
		if "," in af:
			return 0
		af = float(af)
		return(af)
	vcf['AF'] = vcf['INFO'].apply(alt_af)
 
	for eqtl_index in range(len(eqtls)): 

		with open(f"{base_output_dir}/hg19/eqtl/eqtl_peer_sumstats{locus_num}_{settings['peer_fraction_var_explained'][eqtl_index]}.txt", "w") as w:
			w.write("feature\trsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref\talt\teffect_af\tbeta\tse\tzscore\tpvalue\tN\n")
			for i in range(vcf.shape[0]):
				var = vcf.iloc[i, :]
				e = eqtls[eqtl_index][i]
				variant_id = f"{var['#CHROM']}_{var['POS']}_{var['REF']}_{var['ALT']}_hg19"
				w.write(f"nameless_gene\t{var['ID']}\t{variant_id}\t{hg19}\t{var['#CHROM']}\t{var['POS']}\t{var['REF']}\t{var['ALT']}\t{var['AF']}\t{e[0]}\t{e[1]}\t{e[0] / e[1]}\t{e[2]}\t{settings['current_run']['eqtl_sample_size']}\n")

		# bgzip and tabix
		subprocess.run(f'bgzip -f {base_output_dir}/hg19/eqtl/eqtl_peer_sumstats{locus_num}_{settings["peer_fraction_var_explained"][eqtl_index]}.txt'.split())
		subprocess.run(f"tabix -f -S 1 -s 5 -b 6 -e 6 {base_output_dir}/hg19/eqtl/eqtl_peer_sumstats{locus_num}_{settings['peer_fraction_var_explained'][eqtl_index]}.txt.gz".split())


def write_answer_key(gwas_effect_sizes, eqtl_effect_sizes, locus, index):
	with open(f"{base_output_dir}/answer_key.txt", "a") as a:
		# Write GWAS/eQTL variants, effect sizes...might also add more info later
		info = ""
		for i in range(len(gwas_effect_sizes)):
			if gwas_effect_sizes[i] != 0:
				info += "gwas:" + settings["current_run"]["rsids"].iloc[i] + ":" + str(gwas_effect_sizes[i]) + ","
		for i in range(len(eqtl_effect_sizes)):
			if eqtl_effect_sizes[i] != 0:
				info += "eqtl:" + settings["current_run"]["rsids"].iloc[i] + ":" + str(eqtl_effect_sizes[i]) + ","
		a.write(f'{index}\t{locus[0]}\t{locus[1]}\t{info}\t{settings["current_run"]["gwas_case_sample_size"]}\t{settings["current_run"]["gwas_control_sample_size"]}\t{settings["current_run"]["eqtl_sample_size"]}\t{len(gwas_effect_sizes)}\n')

def run_liftover(settings):
	print(f"python scripts/simulate_data/liftover_sumstats_hg19_to_hg38.py {base_output_dir} {tmp_dir} {locus_num}")
	subprocess.run(f"python scripts/simulate_data/liftover_sumstats_hg19_to_hg38.py {base_output_dir} {tmp_dir} {locus_num}".split())

if __name__ == "__main__":
	main()
