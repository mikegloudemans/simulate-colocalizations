# Author: Mike Gloudemans
# Date: 7/23/2018

import glob
import subprocess
import sys
import gzip

chain_file = "/oak/stanford/groups/smontgom/shared/liftOver/chains/hg19ToHg38.over.chain.gz"

base_dir = sys.argv[1]
tmp_dir = sys.argv[2]
locus_num = int(sys.argv[3])

eqtl_file = f"{base_dir}/hg19/eqtl/eqtl_sumstats{locus_num}.txt.gz"
gwas_file = f"{base_dir}/hg19/gwas/gwas_sumstats{locus_num}.txt.gz"

#
# eQTL
#

# Write liftable file
with open(f"{tmp_dir}/liftable.bed", "w") as w:
	with gzip.open(eqtl_file) as f:
		f.readline()
		for line in f:
			data = line.decode('utf-8').strip().split()
			w.write(f"{'chr'+data[4]}\t{data[5]}\t{int(data[5])+1}\t{data[2]}\n")
# Lift it
subprocess.run(f"liftOver {tmp_dir}/liftable.bed {chain_file} {tmp_dir}/lifted.bed {tmp_dir}/unlifted.bed".split())

# Merge back together with the original version to save metadata; scrap unnecessary rows
with open(f"{tmp_dir}/lifted_sorted.bed", "w") as w:
	subprocess.run(f"sort -k4,4 {tmp_dir}/lifted.bed".split(), stdout=w)

ps = subprocess.Popen(f"zcat {eqtl_file}".split(), stdout=subprocess.PIPE)
with open(f"{tmp_dir}/sumstats_sorted.bed", "w") as w:
	subprocess.run("sort -k3,3".split(), stdout=w, stdin = ps.stdout)
ps.wait()
#subprocess.check_call("zcat {0} | sort -k3,3 > /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed".format(eqtl_file), shell=True)
ps1 = subprocess.Popen(f"join -1 4 -2 3 {tmp_dir}/lifted_sorted.bed {tmp_dir}/sumstats_sorted.bed".split(), stdout=subprocess.PIPE)
ps2 = subprocess.Popen(["sed", "s/\ /\\\t/g"], stdin=ps1.stdout, stdout=subprocess.PIPE)
ps3 = subprocess.Popen("sed s/chr//g".split(), stdin=ps2.stdout, stdout=subprocess.PIPE)
with open(f"{tmp_dir}/lifted_and_joined.bed", "w") as w:
	subprocess.run("sort -k3,3n".split(), stdin=ps3.stdout, stdout=w)
ps1.wait()
ps2.wait()
ps3.wait()

#subprocess.check_call("join -1 4 -2 3 /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed | sed s/\ /\\\t/g | sed s/chr//g | sort -k3,3n > /users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed", shell=True)

# Reformat it one more time
with open(f"{base_dir}/hg38/eqtl/eqtl_sumstats{locus_num}.txt", "w") as w:
	w.write("gene\trsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref\talt\teffect_af\teffect_size\tse\tzscore\tpvalue\tN\n")
	with open(f"{tmp_dir}/lifted_and_joined.bed") as f:
		for line in f:
			data = line.strip().split()
			w.write(f"{data[4]}\t{data[5]}\t{data[1]}_{data[2]}_{data[9]}_{data[10]}_hg38\thg38\t{data[1]}\t{data[2]}\t{data[9]}\t{data[10]}\t{data[11]}\t{data[12]}\t{data[13]}\t{data[14]}\t{data[15]}\t{data[16]}\n")

subprocess.run(f"bgzip -f {base_dir}/hg38/eqtl/eqtl_sumstats{locus_num}.txt".split())
subprocess.run(f"tabix -f -S 1 -s 5 -b 6 -e 6 {base_dir}/hg38/eqtl/eqtl_sumstats{locus_num}.txt.gz".split())

# NOTE: Right now, some will fail indexing because chromosomes don't agree. For now
# I'm just going to let that happen and throw these sites away from consideration for
# hg38
			
#
# GWAS
#

# Write liftable file
with open(f"{tmp_dir}/liftable.bed", "w") as w:
	with gzip.open(gwas_file) as f:
		f.readline()
		for line in f:
			data = line.decode('utf-8').strip().split()
			w.write(f"{'chr'+data[3]}\t{data[4]}\t{int(data[4])+1}\t{data[1]}\n")

# Lift it
subprocess.run(f"liftOver {tmp_dir}/liftable.bed {chain_file} {tmp_dir}/lifted.bed {tmp_dir}/unlifted.bed".split())

# Merge back together with the original version to save metadata; scrap unnecessary rows
with open(f"{tmp_dir}/lifted_sorted.bed", "w") as w:
	subprocess.run(f"sort -k4,4 {tmp_dir}/lifted.bed".split(), stdout=w)

ps = subprocess.Popen(f"zcat {gwas_file}".split(), stdout=subprocess.PIPE)
with open(f"{tmp_dir}/sumstats_sorted.bed", "w") as w:
	subprocess.run("sort -k2,2".split(), stdout=w, stdin = ps.stdout)
ps.wait()

ps1 = subprocess.Popen(f"join -1 4 -2 2 {tmp_dir}/lifted_sorted.bed {tmp_dir}/sumstats_sorted.bed".split(), stdout=subprocess.PIPE)
ps2 = subprocess.Popen(["sed", "s/\ /\\\t/g"], stdin=ps1.stdout, stdout=subprocess.PIPE)
ps3 = subprocess.Popen("sed s/chr//g".split(), stdin=ps2.stdout, stdout=subprocess.PIPE)
with open(f"{tmp_dir}/lifted_and_joined.bed", "w") as w:
	subprocess.run("sort -k3,3n".split(), stdin=ps3.stdout, stdout=w)
ps1.wait()
ps2.wait()
ps3.wait()



# Reformat it one more time

with open(f"{base_dir}/hg38/gwas/gwas_sumstats{locus_num}.txt", "w") as w:
	w.write("rsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref\talt\teffect_af\teffect_size\tse\tzscore\tpvalue\tn_cases\tn_controls\n")
	with open(f"{tmp_dir}/lifted_and_joined.bed") as f:
		for line in f:
			data = line.strip().split()
			w.write(f"{data[4]}\t{data[1]}_{data[2]}_{data[8]}_{data[9]}_hg38\thg38\t{data[1]}\t{data[2]}\t{data[8]}\t{data[9]}\t{data[10]}\t{data[11]}\t{data[12]}\t{data[13]}\t{data[14]}\t{data[15]}\t{data[16]}\n")

subprocess.run(f"bgzip -f {base_dir}/hg38/gwas/gwas_sumstats{locus_num}.txt".split())
subprocess.run(f"tabix -f -S 1 -s 4 -b 5 -e 5 {base_dir}/hg38/gwas/gwas_sumstats{locus_num}.txt.gz".split())

# NOTE: Right now, some will fail indexing because chromosomes don't agree. For now
# I'm just going to let that happen and throw these sites away from consideration for
# hg38
