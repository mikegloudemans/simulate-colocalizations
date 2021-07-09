
files = dir("/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/gwas")
files = files[grepl("gwas_sumstats", files)]
files = files[!grepl("tbi", files)]

all_data = list()

for (file in files[1:20])
{
	all_data[[file]] = read.table(paste0("/users/mgloud/projects/coloc_comparisons/output/simulations/2018-07-27_15-23-15/hg19/gwas/", file), header=TRUE)
}

data = do.call(rbind, all_data)

qqplot(-log10(runif(length(data$pvalue), 0, 1)), -log10(data$pvalue), pch=18, ylab="Empirical distribution of GWAS p-values", xlab="Expected distribution of GWAS p-values")
abline(a=0, b=1, col="red", lty=3)

gen_chisq = data$zscore^2
median(gen_chisq)/qchisq(0.5,1)
