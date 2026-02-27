#!/bin/bash
set +o posix
set -eo pipefail

input_vcf=H_ota_filtered_PASS.vcf.gz
sp_name=H_ota
sex_info=sex.deal.txt
cpu=40
## sex_info
## id (tab) sex (1male;2female)
## ZD-1-1   1
## ZD-1-10 2

plink --vcf ${input_vcf} --out ${sp_name}1 --allow-extra-chr --memory 1000 --biallelic-only
plink --bfile ${sp_name}1 --geno 0.02 --mind 0.02 --make-bed --out ${sp_name}2 --allow-extra-chr --allow-no-sex
plink --bfile ${sp_name}2 --maf 0.02 --make-bed --out ${sp_name}3 --allow-extra-chr --allow-no-sex
sed 's/ /\t/g' ${sp_name}3.fam | awk 'NR==FNR{a[$1]=$2}NR!=FNR{print $0"\t"a[$1]}' ${sex_info} - | csvtk -t cut -f 1,2,3,4,7,7 | sed 's/\t/ /g' > ${sp_name}3.fam.deal
mv ${sp_name}3.fam.deal ${sp_name}3.fam

awk '{print $1"\t"$1"_"$4"\t"$3"\t"$4"\t"$5"\t"$6}' ${sp_name}3.bim > ${sp_name}3.bim.deal
mv ${sp_name}3.bim.deal ${sp_name}3.bim

gemma-0.98.5-linux-static-AMD64 -bfile ${sp_name}3 -gk 2 -o kin -r2 0.8
gemma-0.98.5-linux-static-AMD64 -bfile ${sp_name}3 -k output/kin.sXX.txt -lmm 4 -o ${sp_name}_gwas_lmm

ln -s output/${sp_name}_gwas_lmm.assoc.txt

csvtk -t cut -f 1,2,3,15 ${sp_name}_gwas_lmm.assoc.txt | sed '1d' | csvtk cut -t -f 2,2,4 | sed 's/_/\t/2' | csvtk -t add-header -n SNP_id,CHR,BP,P > ${sp_name}_gwas_lmm.assoc.txt.p
python zhuanhuan.py ${sp_name}_gwas_lmm.assoc.txt.p ${sp_name}_gwas_lmm.assoc.txt.p.deal 4

Rscript QQplot.r ${sp_name}_gwas_lmm.assoc.txt.p.deal &> lamda.txt
## change significant value here
Rscript Manhattan.r ${sp_name}_gwas_lmm.assoc.txt.p.deal gwas_lmm.mhd 4.35E-5 2.17E-3

rm ${sp_name}1* ${sp_name}2* ${sp_name}_gwas_lmm.assoc.txt.p

awk '$2=="1"' ${sex_info} | awk '{print $1}' > male.list
awk '$2=="2"' ${sex_info} | awk '{print $1}' > female.list

bcftools view --threads ${cpu} -S male.list ${input_vcf} -O z -o males.vcf.gz
bcftools view --threads ${cpu} -S female.list ${input_vcf} -O z -o females.vcf.gz

vcftools --gzvcf males.vcf.gz --hardy --out males_hardy
vcftools --gzvcf females.vcf.gz --hardy --out females_hardy

python calc_sex_het.py
python find_sdr_region.py
sort -k 5nr,5 sex_heterozygosity_diff.csv > sex_heterozygosity_diff.sort.tsv

rm males.vcf.gz females.vcf.gz sex_heterozygosity_diff.csv

pigz --best ${sp_name}_gwas_lmm.assoc.txt.p.deal sex_heterozygosity_diff.sort.csv
