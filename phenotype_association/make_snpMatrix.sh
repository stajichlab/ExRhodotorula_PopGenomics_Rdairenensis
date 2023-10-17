#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 8gb 
module load pyvcf
VCF=Rdar_v1.All.polymorphicDBVPG_filtered.vcf.gz
REFGENOME=../genome/Rhodotorula_dairenensis_NRRL_Y-2504.scaffolds.fa
OUTMATRIX=$(basename $VCF .vcf.gz).matrix.tsv
../scripts/snpEff_2_tab.py $VCF $REFGENOME > $OUTMATRIX
