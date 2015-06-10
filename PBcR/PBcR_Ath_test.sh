#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l mem=64gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR


start=`date +%s`

module load java/1.8.0_25

PBcR_bin=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/wgs-8.3rc2/Linux-amd64/bin
spec_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Arabidopsis/pacbio.test.spec
fq_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Arabidopsis/Pacbio.fastq
lib_name=Ath
gsize=135000000

#pacbio.spec
## limit to 32GB. By default the pipeline will auto-detect memory and try to use maximum. This allow limiting it
#ovlMemory = 32
#ovlStoreMemory= 32000
#merylMemory = 32000

#fast consensus (-pbCNS parameter)
$PBcR_bin/PBcR -threads $PBS_NP -length 500 -partitions 200 -l $lib_name -s $spec_file -fastq $fq_file genomeSize=$gsize

#ref=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Lambda_Phase/sampleData/reference.fasta
#asm=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/PBcR/lambda/9-terminator/asm.ctg.fasta
#/usr/local/bin/dnadiff $ref $asm

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

