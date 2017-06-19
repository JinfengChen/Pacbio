#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=stdout
#SBATCH -p intel
#SBATCH --workdir=./

GATK=/opt/linux/centos/7.x/x86_64/pkgs/gatk/3.4-46/GenomeAnalysisTK.jar
JAVA=/opt/linux/centos/7.x/x86_64/pkgs/java/jdk1.7.0_17/bin/java
genome=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Fairchild/Fairchild_v1.fasta
CPU=$SLURM_NTASKS

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

#echo "ReduceReads for each bam using nway co-reduce, which could keep these reads that not have SNP in some sample in common SNPs site"
#$JAVA -Xmx10g -jar $GATK -T ReduceReads -R $genome -I $BAM -o reduced.bam
#echo "Reduce Done"
echo "Single-sample call SNP from illumina and Pacbio reads"

#-L scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8,scaffold_9
#files=(LAP CHP SSO SWO PKM CLM WLM WMM FAIRCHILD)
#files=(SWO.30X)
bam=`ls *.bam | head -n $N | tail -n 1`
prefix=${bam%.bam}
#snp call
if [ ! -e $prefix.gatk.raw.g.vcf ]; then
$JAVA -Xmx40g -jar $GATK \
      -T HaplotypeCaller \
      -R $genome \
      -I $bam \
      -o $prefix.gatk.raw.g.vcf \
      -nct $CPU \
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
