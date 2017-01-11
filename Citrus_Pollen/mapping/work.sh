echo "run mapping"
qsub -t 1-7 run_speedseq_qsub.sh
mv *.discordants.* *.splitters.* trash/

echo "summary bam"
python Run_Qualimap.py --bam ./
python Sum_Qualimap.py --bam ./ > Citrus.bam.summary

echo "Pollen SNP"

grep "^scaffold" Pollen41.gatk.snp.hardfilter.vcf| grep "PASS" | awk '$10~/1\/1/' | wc -l
grep "^scaffold" Pollen41.gatk.snp.hardfilter.vcf| grep "PASS" | awk '$10~/0\/1/' | wc -l

echo "genotype"
cat /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Citrus_Pollen/fairchild_SNP/Fairchild.vcf.recode.minDP5.NonRef.recode.vcf | awk 'BEGIN{FS="\t"}($0 !~ /^#/){if(length($5)<= 1 && length($4)<= 1) print $0} ($0 ~/^#/){print $0}' > dbSNP.clean.vcf
python removePGT_PID.py --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Citrus_Pollen/fairchild_SNP/Fairchild.vcf.recode.minDP5.NonRef.recode.vcf > dbSNP.clean.vcf &

