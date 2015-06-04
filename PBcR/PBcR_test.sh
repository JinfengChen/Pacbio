#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l mem=4gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR


start=`date +%s`

module load java/1.8.0_25

PBcR_bin=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/wgs-8.3rc2/Linux-amd64/bin
spec_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Lambda_Phase/sampleData/pacbio.spec
fq_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Lambda_Phase/sampleData/pacbio.filtered_subreads.fastq
lib_name=lambda
gsize=50000

$PBcR_bin/PBcR -threads 2 -length 500 -partitions 200 -l $lib_name -s $spec_file -fastq $fq_file genomeSize=$gsize

ref=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Lambda_Phase/sampleData/reference.fasta
asm=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/PBcR/lambda/9-terminator/asm.ctg.fasta
#/usr/local/bin/dnadiff $ref $asm

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

