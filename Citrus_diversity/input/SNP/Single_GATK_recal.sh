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
BASE=LAP

if [ ! -f $BASE.recal_data.grp ] ; then 
 echo "Generating $BASE.recal_data.grp from $BASE.bam"
 java -Xmx15g -jar $GATK \
 -T BaseRecalibrator \
 -nct $cpu \
 -I $BASE.bam \
 -R $genome \
 -knownSites LAP.freebayes.vcf \
 -o $BASE.recal_data.grp 
fi

if [ ! -f $BASE.recal.bam ] ; then 
  echo "Generating $BASE.recal.bam from $BASE.bam"
  java -Xmx15g -jar $GATK \
   -T PrintReads \
   -nct $cpu \
   -R $genome \
   -I $BASE.bam \
   -BQSR $BASE.recal_data.grp \
   -o $BASE.recal.bam
fi


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

