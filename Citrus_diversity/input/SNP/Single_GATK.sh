#!/bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=50gb
#PBS -l walltime=200:00:00
#PBS -d ./

#module load gatk/3.4-46
#GATK=/opt/GATK/2.4-3-g2a7af43/GenomeAnalysisTK.jar
GATK=/opt/linux/centos/7.x/x86_64/pkgs/gatk/3.4-46/GenomeAnalysisTK.jar
Picard=/opt/picard/1.81/
JAVA=/opt/linux/centos/7.x/x86_64/pkgs/java/jdk1.7.0_17/bin/java
samtools=/usr/bin/samtools
genome=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/SNP_calling/input/Cclementina_v1.0_scaffolds.fa
#dbSNP=/rhome/cjinfeng/HEG4_cjinfeng/Variations/dbSNP_rice/dbSNP_HEG4/HEG4_dbSNP.vcf
repeat=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/SNP_calling/input/Clementina.fasta.RepeatMasker.out.bed
BAM=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/SNP_calling/input/Citrus_asm40_preads.bam
prefix=PKM
bam=$prefix\.bam
cpu=$PBS_NP

#echo "ReduceReads for each bam using nway co-reduce, which could keep these reads that not have SNP in some sample in common SNPs site"
#$JAVA -Xmx10g -jar $GATK -T ReduceReads -R $genome -I $BAM -o reduced.bam
#echo "Reduce Done"
echo "Single-sample call SNP from illumina and Pacbio reads"

:<<SKIP
###Markduplicate
if [ ! -e $prefix.dedup.bam ]; then
$JAVA -Xmx10g -jar $Picard/MarkDuplicates.jar \
      ASSUME_SORTED=TRUE \
      REMOVE_DUPLICATES=TRUE \
      VALIDATION_STRINGENCY=LENIENT \
      INPUT=$BAM \
      OUTPUT=$prefix.dedup.bam \
      METRICS_FILE=$prefix.dedup.metrics
fi

###Read Group
if [ ! -e $prefix.rg.bam ]; then
$JAVA -Xmx10g -jar $Picard/AddOrReplaceReadGroups.jar \
      INPUT=$prefix.dedup.bam \
      OUTPUT=$prefix.rg.bam \
      SORT_ORDER=coordinate \
      CREATE_INDEX=true \
      RGID=Fairchild \
      RGLB=smt40 \
      RGPU=AGI \
      RGPL=Pacbio \
      RGSM=Citrus \
      VALIDATION_STRINGENCY=SILENT

$samtools flagstat $prefix.rg.bam > $prefix.rg.bam.flagstat
fi
SKIP

#snp call
if [ ! -e $prefix.gatk.snp.raw.gvcf ]; then
$JAVA -Xmx40g -jar $GATK \
      -T HaplotypeCaller \
      -R $genome \
      -I $bam \
      -o $prefix.gatk.raw.g.vcf \
      -nct $PBS_NP \
      -allowPotentiallyMisencodedQuals \
      --allow_potentially_misencoded_quality_scores \
      --genotyping_mode DISCOVERY \
      --emitRefConfidence GVCF \
      -stand_call_conf 30 \
      -stand_emit_conf 10 \
      --defaultBaseQualities 30 \
      -variant_index_type LINEAR \
      -variant_index_parameter 128000
fi

#$JAVA -Xmx10g -jar $GATK \
#      -T SelectVariants \
#      -R $genome \
#      -V $prefix.gatk.raw.vcf \
#      -selectType SNP \
#      -o $prefix.gatk.snp.raw.vcf

#$JAVA -Xmx10g -jar $GATK \
#      -T SelectVariants \
#      -R $genome \
#      -V $prefix.gatk.raw.vcf \
#      -selectType INDEL \
#      -o $prefix.gatk.indel.raw.vcf

fi

:<<SKIP
###hardfilter indel
if [ ! -e $prefix.gatk.indel.pass.vcf ]; then
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.indel.raw.vcf \
      -o $prefix.gatk.indel.hardfilter.vcf \
      --filterExpression "QD < 2.0" \
      --filterName "QDFilter" \
      --filterExpression "ReadPosRankSum < -20.0" \
      --filterName "ReadPosFilter" \
      --filterExpression "FS > 200.0" \
      --filterName "FSFilter" \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --filterExpression "QUAL < 30.0 || DP < 6 || DP > 5000 || HRun > 5" \
      --filterName "QualFilter"

$JAVA -Xmx1g -jar $GATK -T SelectVariants -R $genome --variant $prefix.gatk.indel.hardfilter.vcf -o $prefix.gatk.indel.pass.vcf --excludeFiltered

fi

###hardfilter snp
if [ ! -e $prefix.gatk.snp.pass.vcf ]; then
$JAVA -Xmx1g -jar $GATK \
      -T VariantFiltration \
      -R $genome \
      --variant $prefix.gatk.snp.raw.vcf \
      -o $prefix.gatk.snp.hardfilter.vcf \
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
      --mask $prefix.gatk.indel.pass.vcf \
      --maskName "INDEL"

$JAVA -Xmx1g -jar $GATK -T SelectVariants -R $genome --variant $prefix.gatk.snp.hardfilter.vcf -o $prefix.gatk.snp.pass.vcf --excludeFiltered
fi

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

echo "Done!"
