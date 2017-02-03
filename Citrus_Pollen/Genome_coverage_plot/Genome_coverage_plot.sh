#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=30G
#SBATCH --time=10:00:00
#SBATCH --output=Genome_coverage_plot.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#sbatch --array 1-2 Genome_coverage_plot.sh

module load samtools/0.1.19
#module load bedtools/2.25.0
bedtools=~/BigData/software/bedtools2-2.19.0/bin/bedtools
PATH=$PATH:~/BigData/software/SVcaller/ROOT/bin/
chrlen=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.9.chrlen
window=100000
step=100000

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

FILE=`ls Pollen41.bam | head -n $N | tail -n 1`
SAMPLE=${FILE%.bam}
echo "File: $FILE"
echo "Sample: $SAMPLE"
prefix=$SAMPLE

if [ ! -e $SAMPLE\.filter_MQ2.bam ]; then
    echo "Filtering $SAMPLE ..."
    samtools view -@ $CPU -bq 2 $SAMPLE\.bam > $SAMPLE\.filter_MQ2.bam
fi
$bedtools makewindows -g $chrlen -w $window -s $step > $prefix\.genome.slidingwin_$window\_$step\.bed
$bedtools coverage -abam $SAMPLE\.filter_MQ2.bam -b $prefix\.genome.slidingwin_$window\_$step\.bed -hist > $prefix\.genome.slidingwin_$window\_$step\.bed.hist.cov
python Bam_Cov_SlidingWin.py --input $prefix\.genome.slidingwin_$window\_$step\.bed.hist.cov --window $window --step $step
python Convert_position_to_one_chromosome_bed.py --chr $chrlen --bed $prefix\.genome.slidingwin_$window\_$step\.bed.hist.cov.bed

echo "Done"
