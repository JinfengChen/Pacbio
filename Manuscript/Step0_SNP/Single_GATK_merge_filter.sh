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
prefix=Citrus_9strain

if [ ! -f $merge_vcf ] ; then
$JAVA -Xmx15g -jar $GATK \
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

####split SNP and indel
if [ ! -e $prefix.gatk.INDEL.raw.vcf ]; then
$JAVA -Xmx10g -jar $GATK \
   -R $genome \
   -T SelectVariants \
   -V $merge_vcf \
   -o $prefix.gatk.INDEL.raw.vcf \
   -selectType INDEL
fi

if [ ! -e $prefix.gatk.SNP.raw.vcf ]; then
$JAVA -Xmx10g -jar $GATK \
   -R $genome \
   -T SelectVariants \
   -V $merge_vcf \
   -o $prefix.gatk.SNP.raw.vcf \
   -selectType SNP
fi


###hardfilter indel
if [ ! -e $prefix.gatk.INDEL.pass.vcf ]; then
#$JAVA -Xmx10g -jar $GATK \
#      -T VariantFiltration \
#      -R $genome \
#      --variant $prefix.gatk.INDEL.raw.vcf \
#      -o $prefix.gatk.INDEL.hardfilter.vcf \
#      --filterExpression "QD < 2.0" \
#      --filterName "QDFilter" \
#      --filterExpression "ReadPosRankSum < -20.0" \
#      --filterName "ReadPosFilter" \
#      --filterExpression "FS > 200.0" \
#      --filterName "FSFilter" \
#      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
#      --filterName "HARD_TO_VALIDATE" \
#      --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
#      --filterName "QualFilter"

$JAVA -Xmx1g -jar $GATK -T SelectVariants -R $genome --variant $prefix.gatk.INDEL.hardfilter.vcf -o $prefix.gatk.INDEL.pass.vcf --excludeFiltered

fi


###hardfilter snp
if [ ! -e $prefix.gatk.SNP.pass.vcf ]; then
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.SNP.raw.vcf \
      -o $prefix.gatk.SNP.hardfilter.vcf \
      --clusterSize 3 \
      --clusterWindowSize 10 \
      --filterExpression "QD < 2.0" \
      --filterName "QDFilter" \
      --filterExpression "MQ < 40.0" \
      --filterName "MQFilter" \
      --filterExpression "FS > 60.0" \
      --filterName "FSFilter" \
      --filterExpression "HaplotypeScore > 13.0" \
      --filterName "HaplotypeScoreFilter" \
      --filterExpression "MQRankSum < -12.5" \
      --filterName "MQRankSumFilter" \
      --filterExpression "ReadPosRankSum < -8.0" \
      --filterName "ReadPosRankSumFilter" \
      --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
      --filterName "StandardFilters" \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --mask $prefix.gatk.INDEL.pass.vcf \
      --maskName "INDEL"
      --mask $repeat
      --maskName "REPEAT"

$JAVA -Xmx1g -jar $GATK -T SelectVariants -R $genome --variant $prefix.gatk.SNP.hardfilter.vcf -o $prefix.gatk.SNP.pass.vcf --excludeFiltered
fi

:<<SKIP
###mask SNPs cluster and INDEL
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.snp.VQSR.vcf \
      -o $prefix.gatk.snp.VQSR.INDEL.vcf \
      --clusterSize 3 \
      --clusterWindowSize 10 \
      --mask $prefix.gatk.indel.pass.vcf \
      --maskExtension 10 \
      --maskName "INDEL"

###mask repeat
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.snp.VQSR.INDEL.vcf \
      -o $prefix.gatk.snp.VQSR.MASK.vcf \
      --mask $repeat \
      --maskName "REPEAT"

$JAVA -Xmx1g -jar $GATK -T SelectVariants -R $genome --variant $prefix.gatk.snp.VQSR.MASK.vcf -o $prefix.gatk.snp.VQSR.pass.vcf --excludeFiltered
SKIP

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

