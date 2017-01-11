#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=10:00:00
#PBS -d ./
#PBS -j oe

#cd $PBS_O_WORKDIR
#qsub -t 1-7 indexbam_qsub.sh

module load samtools

CPU=$PBS_NP
if [ ! $CPU ]; then
   CPU=2
fi

N=$PBS_ARRAYID
if [ ! $N ]; then
    N=$1
fi

FILE=`ls -1 *.bam | head -n $N | tail -n 1`

prefix=${FILE%.*}
if [ ! -e $prefix\.flagstats ]; then
   samtools flagstat $FILE > $prefix\.flagstat
fi


echo "Done"
