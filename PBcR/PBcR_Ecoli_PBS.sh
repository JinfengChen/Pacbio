#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=36gb
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR
##PBS dependency job have problem when doing assembly. do not use that for assembly.

start=`date +%s`

module load java/1.8.0_25

PBcR_bin=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/wgs-8.3rc2/Linux-amd64/bin
spec_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Ecoli/selfSampleData/pacbio.PBS.spec
fq_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Ecoli/selfSampleData/pacbio_filtered.fastq
lib_name=K12
gsize=4650000

#pacbio.spec
## limit to 32GB. By default the pipeline will auto-detect memory and try to use maximum. This allow limiting it
#ovlMemory = 32
#ovlStoreMemory= 32000
#merylMemory = 32000

#PBS_NP=4
#fast consensus (-pbCNS parameter)
$PBcR_bin/PBcR -threads $PBS_NP -length 500 -partitions 200 -l $lib_name -s $spec_file -fastq $fq_file genomeSize=$gsize

#ref=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Ecoli/selfSampleData/reference.fasta
#asm=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/PBcR/K12/9-terminator/asm.ctg.fasta
#/usr/local/bin/dnadiff $ref $asm

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

