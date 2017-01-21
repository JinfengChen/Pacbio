#Fairchild
grep "^scaffold" Fairchild.vcf.recode.minDP5.NonRef.recode.vcf | grep "PASS" | awk '$10~/1\/1/' | wc -l
#401092
grep "^scaffold" Fairchild.vcf.recode.minDP5.NonRef.recode.vcf | grep "PASS" | awk '$10~/0\/1/' | wc -l
#888246
grep "^scaffold" Fairchild.vcf.recode.minDP5.NonRef.recode.vcf | grep "PASS" -c
#1305204

#Pollen41
grep "^scaffold" Pollen41.gatk.snp.pass.vcf | grep "PASS" | awk '$10~/1\/1/' | wc -l
#1089389
grep "^scaffold" Pollen41.gatk.snp.pass.vcf | grep "PASS" | awk '$10~/0\/1/' | wc -l
#220427
grep "^scaffold" Pollen41.gatk.snp.pass.vcf | grep "PASS" -c
#1311019

vcftools --vcf Fairchild.vcf.recode.minDP5.Het.recode.vcf --diff Pollen41.gatk.snp.pass.vcf --diff-site --out Fairchild_v_pollen41 --chr scaffold_1 --chr scaffold_2 --chr scaffold_3 --chr scaffold_4 --chr scaffold_5 --chr scaffold_6 --chr scaffold_7 --chr scaffold_8 --chr scaffold_9
vcftools --vcf Fairchild.vcf.recode.minDP5.Het.recode.vcf --diff Pollen41.gatk.snp.pass.vcf --diff-site --out FairchildHet_v_pollen41 --chr scaffold_1 --chr scaffold_2 --chr scaffold_3 --chr scaffold_4 --chr scaffold_5 --chr scaffold_6 --chr scaffold_7 --chr scaffold_8 --chr scaffold_9
vcftools --vcf Fairchild.vcf.recode.minDP5.Het.recode.vcf --diff Fairchild.gatk.snp.pass.vcf --diff-site --out Fairchild_v_FCillumina --chr scaffold_1 --chr scaffold_2 --chr scaffold_3 --chr scaffold_4 --chr scaffold_5 --chr scaffold_6 --chr scaffold_7 --chr scaffold_8 --chr scaffold_9 > log 2>&1 &


echo "genotype"
awk -F'\t' '(substr($10,1,3)=="0/1")' < dbSNP.clean.vcf > dbSNP.clean_Het.vcf
vcftools --vcf dbSNP.clean_Het.vcf --diff Pollen41_GT.gatk.snp.pass.vcf --diff-site --out FairchildHet_v_pollen41_GT --chr scaffold_1 --chr scaffold_2 --chr scaffold_3 --chr scaffold_4 --chr scaffold_5 --chr scaffold_6 --chr scaffold_7 --chr scaffold_8 --chr scaffold_9
vcftools --vcf dbSNP.clean_Het.vcf --diff Pollen43_GT.gatk.snp.pass.vcf --diff-site --out FairchildHet_v_pollen43_GT --chr scaffold_1 --chr scaffold_2 --chr scaffold_3 --chr scaffold_4 --chr scaffold_5 --chr scaffold_6 --chr scaffold_7 --chr scaffold_8 --chr scaffold_9
vcftools --vcf Pollen43_GT.gatk.snp.pass.vcf --diff Pollen41_GT.gatk.snp.pass.vcf --diff-site --out Pollen43_v_pollen41_GT --chr scaffold_1 --chr scaffold_2 --chr scaffold_3 --chr scaffold_4 --chr scaffold_5 --chr scaffold_6 --chr scaffold_7 --chr scaffold_8 --chr scaffold_9

