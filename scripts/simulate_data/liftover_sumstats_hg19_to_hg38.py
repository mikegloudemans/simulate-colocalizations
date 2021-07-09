# Author: Mike Gloudemans
# Date: 7/23/2018

import glob
import subprocess
import sys
import gzip

base_dir = sys.argv[1]
total_test_sites = int(sys.argv[2])

for i in range(total_test_sites):
    print i
    eqtl_file = "{0}/hg19/eqtl/eqtl_sumstats{1}.txt.gz".format(base_dir, i)
    gwas_file = "{0}/hg19/gwas/gwas_sumstats{1}.txt.gz".format(base_dir, i)

    #
    # eQTL
    #

    # Write liftable file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/liftable.bed", "w") as w:
        with gzip.open(eqtl_file) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                w.write("{0}\t{1}\t{2}\t{3}\n".format("chr"+data[4], data[5], int(data[5])+1, data[2]))
    # Lift it
    subprocess.check_call("liftOver /users/mgloud/projects/coloc_comparisons/tmp/liftable.bed /users/mgloud/software/liftOver/chains/hg19ToHg38.over.chain.gz /users/mgloud/projects/coloc_comparisons/tmp/lifted.bed /users/mgloud/projects/coloc_comparisons/tmp/unlifted.bed", shell=True)
    
    # Merge back together with the original version to save metadata; scrap unnecessary rows
    subprocess.check_call("sort -k4,4 /users/mgloud/projects/coloc_comparisons/tmp/lifted.bed > /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed", shell=True)
    subprocess.check_call("zcat {0} | sort -k3,3 > /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed".format(eqtl_file), shell=True)
    subprocess.check_call("join -1 4 -2 3 /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed | sed s/\ /\\\t/g | sed s/chr//g | sort -k3,3n > /users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed", shell=True)

    # Reformat it one more time
    with open("{0}/hg38/eqtl/eqtl_sumstats{1}.txt".format(base_dir, i), "w") as w:
        w.write("gene\trsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref_allele\talt_allele\teffect_af\teffect_size\tse\tzscore\tpvalue\tN\n")
        with open("/users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed") as f:
            for line in f:
                data = line.strip().split()
                w.write("{11}\t{0}\t{1}_{2}_{3}_{4}_hg38\thg38\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(data[5], data[1], data[2], data[9], data[10], data[11], data[12], data[13], data[14], data[15], data[16], data[4]))

    subprocess.call("bgzip -f {0}/hg38/eqtl/eqtl_sumstats{1}.txt".format(base_dir, i), shell=True)
    subprocess.call("tabix -f -S 1 -s 5 -b 6 -e 6 {0}/hg38/eqtl/eqtl_sumstats{1}.txt.gz".format(base_dir, i), shell=True)

    # NOTE: Right now, some will fail indexing because chromosomes don't agree. For now
    # I'm just going to let that happen and throw these sites away from consideration for
    # hg38
                
    #
    # GWAS
    #

    # Write liftable file
    with open("/users/mgloud/projects/coloc_comparisons/tmp/liftable.bed", "w") as w:
        with gzip.open(gwas_file) as f:
            f.readline()
            for line in f:
                data = line.strip().split()
                w.write("{0}\t{1}\t{2}\t{3}\n".format("chr"+data[3], data[4], int(data[4])+1, data[1]))
    # Lift it
    subprocess.check_call("liftOver /users/mgloud/projects/coloc_comparisons/tmp/liftable.bed /users/mgloud/software/liftOver/chains/hg19ToHg38.over.chain.gz /users/mgloud/projects/coloc_comparisons/tmp/lifted.bed /users/mgloud/projects/coloc_comparisons/tmp/unlifted.bed", shell=True)
    
    # Merge back together with the original version to save metadata; scrap unnecessary rows
    subprocess.check_call("sort -k4,4 /users/mgloud/projects/coloc_comparisons/tmp/lifted.bed > /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed", shell=True)
    subprocess.check_call("zcat {0} | sort -k2,2 > /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed".format(gwas_file), shell=True)
    subprocess.check_call("join -1 4 -2 2 /users/mgloud/projects/coloc_comparisons/tmp/lifted_sorted.bed /users/mgloud/projects/coloc_comparisons/tmp/sumstats_sorted.bed | sed s/\ /\\\t/g  | sort -k3,3n > /users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed", shell=True)

    # Reformat it one more time

    with open("{0}/hg38/gwas/gwas_sumstats{1}.txt".format(base_dir, i), "w") as w:
        w.write("rsid\tvariant_id\tgenome_build\tchr\tsnp_pos\tref_allele\talt_allele\teffect_af\teffect_size\tse\tzscore\tpvalue\tn_cases\tn_controls\n")
        with open("/users/mgloud/projects/coloc_comparisons/tmp/lifted_and_joined.bed") as f:
            for line in f:
                data = line.strip().split()
                w.write("{0}\t{1}_{2}_{3}_{4}_hg38\thg38\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(data[4], data[1], data[2], data[8], data[9], data[10], data[11], data[12], data[13], data[14], data[15], data[16]))

    subprocess.call("bgzip -f {0}/hg38/gwas/gwas_sumstats{1}.txt".format(base_dir, i), shell=True)
    subprocess.call("tabix -f -S 1 -s 4 -b 5 -e 5 {0}/hg38/gwas/gwas_sumstats{1}.txt.gz".format(base_dir, i), shell=True)

    # NOTE: Right now, some will fail indexing because chromosomes don't agree. For now
    # I'm just going to let that happen and throw these sites away from consideration for
    # hg38
