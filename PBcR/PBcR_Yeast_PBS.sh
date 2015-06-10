#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=64gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR


start=`date +%s`

module load java/1.8.0_25

PBcR_bin=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/wgs-8.3rc2/Linux-amd64/bin
spec_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/pacbio.PBS.spec
fq_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/Pacbio.fastq
lib_name=W303
gsize=12495682

#pacbio.spec
## limit to 32GB. By default the pipeline will auto-detect memory and try to use maximum. This allow limiting it
#ovlMemory = 32
#ovlStoreMemory= 32000
#merylMemory = 32000

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

