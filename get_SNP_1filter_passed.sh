#!/bin/bash

## Do Fisher, Chi2 and Inbred Coefficient
Rscript fisher_quality_test.R Count_Geno-$1-all.vcf.txt
Rscript chi2_quality_test.R Count_Geno-$1-all.vcf.txt
Rscript inbred_quality_test.R stat-$1.txt

# remove ".", NA and >2 allele SNPs from inbred
awk 'BEGIN{OFS="\t"}{print $1"_"$2,$1,$2,$3}' inbred_quality_stat-$1.txt > modified_inbred_quality_stat-$1.txt
grep -v "^#" $1-all.vcf  | awk '$4=="." || length($4)>1 || $5=="." || length($5)>1 {print $3}' | cut -c 2- > remove_these_pos.txt
grep -Fwv -f remove_these_pos.txt modified_inbred_quality_stat-$1.txt | awk 'BEGIN{OFS="\t"}{print $2,$3,$4}' > inbred_quality_stat-$1.txt
rm remove_these_pos.txt modified_inbred_quality_stat-$1.txt

## passed locations
awk -F"\t" '$3 > 3' fisher_quality_Count_Geno-$1-all.vcf.txt > passed_fisher.txt
awk -F"\t" '$3 > 90' chi2_quality_Count_Geno-$1-all.vcf.txt > passed_chi2.txt
awk -F"\t" '$3 > 80' inbred_quality_stat-$1.txt > passed_inbred.txt

## remove header
sed -i 1d passed_fisher.txt
sed -i 1d passed_chi2.txt
sed -i 1d passed_inbred.txt

# get location like "1A_100000" to grep all.vcf file
awk -F"\t" '{print $1"_"$2}' passed_fisher.txt > pos_passed_fisher.txt
awk -F"\t" '{print $1"_"$2}' passed_chi2.txt > pos_passed_chi2.txt
awk -F"\t" '{print $1"_"$2}' passed_inbred.txt > pos_passed_inbred.txt

# get vcf files passing each filter
grep -m 1 "^#CHR" $1-all.vcf > $1.Fisher.vcf
grep -m 1 "^#CHR" $1-all.vcf > $1.chi2.vcf
grep -m 1 "^#CHR" $1-all.vcf > $1.inbred.vcf
grep -F -f pos_passed_fisher.txt  $1-all.vcf >> $1.Fisher.vcf
grep -F -f pos_passed_chi2.txt  $1-all.vcf >> $1.chi2.vcf
grep -F -f pos_passed_inbred.txt  $1-all.vcf >> $1.inbred.vcf

# Recover SNP positions passing at least one filtering criteria:
mkdir 01_passed_SNPs
mv *vcf 01_passed_SNPs
cd 01_passed_SNPs
/homes/sshrest1/scripts/tassel/extract_passed_SNPs.sh $1

mkdir 01_snps_1filter_passed 02_other_vcf_files 03_SNP_info 
mv *SNPs-1filter-passed.vcf 01_snps_1filter_passed
mv *vcf 02_other_vcf_files
mv all_pos_from_three_filters.txt FILTER_info.txt filter_passed.txt 03_SNP_info

cd ..
#mkdir 03_tmp_files 
#mv passed_fisher.txt passed_chi2.txt passed_inbred.txt pos_passed_* 03_tmp_files
rm passed_fisher.txt passed_chi2.txt passed_inbred.txt pos_passed_*

mkdir 02_quality_files
mv chi* inbred* stat* fisher* Count_Geno-$1-all.vcf.txt 02_quality_files






