#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l mem=100gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./


#cd $PBS_O_WORKDIR
#qsub -t 1-7 indexbam_qsub.sh

module load samtools
PATH=$PATH:~/BigData/software/SVcaller/ROOT/bin/
genome=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa

CPU=$PBS_NP
if [ ! $CPU ]; then
   CPU=2
fi

N=$PBS_ARRAYID
if [ ! $N ]; then
    N=$1
fi

FILE=`ls *_1.fastq.gz | grep _1\.fastq\.gz | head -n $N | tail -n 1`
R1=$FILE
R2=`echo $R1 | perl -p -e 's/_1\.fastq/_2.fastq/'`
SAMPLE=${FILE%_1.fastq.gz}

if [ ! -e $SAMPLE\.bam ]; then
echo "mapping $SAMPLE ..."
/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq_20161206/speedseq/bin/speedseq align \
     -t $PBS_NP \
     -o $SAMPLE \
     -R "@RG\tID:id\tSM:$SAMPLE\tLB:lib" \
     $genome \
     $SAMPLE\_1.fastq.gz \
     $SAMPLE\_2.fastq.gz
fi

echo "Done"
