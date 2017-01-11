#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=20g
#PBS -l walltime=100:00:00
#PBS -d ./
#PBS -j oe

#cd $PBS_O_WORKDIR

module load bamtools
module load samtools

prefix=Pollen43
p=0.86
bam=$prefix\.bam
bam_out=$prefix\.60X.bam
Picard=/opt/picard/1.81/

/opt/linux/centos/7.x/x86_64/pkgs/picard/2.6.0/bin/picard DownsampleSam I=$bam O=$bam_out P=$p
#samtools sort A119.HEG4allpathv1_BWA.ALL.random.bam A119.HEG4allpathv1_BWA.ALL.random.sort


echo "Done"
