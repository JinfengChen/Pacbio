#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

GATK=/opt/linux/centos/7.x/x86_64/pkgs/gatk/3.4-46/GenomeAnalysisTK.jar
Picard=/opt/picard/1.81/
JAVA=/opt/linux/centos/7.x/x86_64/pkgs/java/jdk1.7.0_17/bin/java
samtools=/usr/bin/samtools
genome=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/SNP_calling/input/Cclementina_v1.0_scaffolds.fa
#dbSNP=/rhome/cjinfeng/HEG4_cjinfeng/Variations/dbSNP_rice/dbSNP_HEG4/HEG4_dbSNP.vcf
repeat=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/SNP_calling/input/Clementina.fasta.RepeatMasker.out.bed
BAM=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/SNP_calling/input/Citrus_asm40_preads.bam
cpu=$PBS_NP
merge_vcf=Citrus_9strain.vcf

if [ ! -f $merge_vcf ] ; then
  java -Xmx15g -jar $GATK \
   -T GenotypeGVCFs \
   -nt $cpu \
   -V CHP.gatk.raw.g.vcf \
   -V LAP.gatk.raw.g.vcf \
   -V WMM.gatk.raw.g.vcf \
   -V WLM.30X.gatk.raw.g.vcf \
   -V CLM.30X.gatk.raw.g.vcf \
   -V PKM.gatk.raw.g.vcf \
   -V SSO.gatk.raw.g.vcf \
   -V SWO.30X.gatk.raw.g.vcf \
   -V FAIRCHILD.gatk.raw.g.vcf \
   -R $genome \
   -o $merge_vcf
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

