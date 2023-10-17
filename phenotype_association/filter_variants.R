#!/usr/bin/env Rscript

library(vcfR)

vcf <- read.vcfR("phenotype_association/Rdar_v1.All.snpEff.vcf.gz")
dna <- ape::read.dna("phenotype_association/Rhodotorula_dairenensis_NRRL_Y-2504.scaffolds.fa.gz", format = "fasta")
gff <- read.table("phenotype_association/Rhodotorula_dairenensis_NRRL_Y-2504.gff3.gz", sep="\t", quote="")
chrom <- create.chromR(name='scaffold', vcf=vcf, seq=dna, ann=gff)
plot(chrom)
chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59.9,  max_MQ = 60.1)
plot(chrom)
chrom <- proc.chromR(chrom, verbose=TRUE)
plot(chrom)
chromoqc(chrom, dp.alpha=20)

head(is.polymorphic(vcf, na.omit = TRUE))
dbvpgonly <- vcf[,1:8]

ad <- extract.gt(dbvpgonly, element = 'AD')
knitr::kable(ad[c(1:2,11,30),1:7])
polyMorphicsites <- is.polymorphic(dbvpgonly, na.omit = TRUE)
polyMorphicDBVPG <- dbvpgonly[polyMorphicsites,]

write.vcf(polyMorphicDBVPG, "phenotype_association/Rdar_v1.All.polymorphicDBVPG_filtered.vcf.gz")
