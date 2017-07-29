#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=run_speedseq_qsub.sh.%A_%a.stdout
#SBATCH -p intel
#SBATCH --workdir=./

#sbatch --array 1 run_speedseq_qsub.sh

#module unload perl/5.20.2
#module load perl/5.14.2
#module load PASA/2.0.2
#module load gmap/2014.12.28
#module load kent/318
#export PERL5LIB=/opt/linux/centos/7.x/x86_64/pkgs/perl/5.20.2/lib/site_perl/5.20.2/x86_64-linux/:$PERL5LIB

start=`date +%s`

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
cpu=12 # need to between 1 and 16
#seqclean transcripts.fasta -c $cpu
#PASA=/opt/linux/centos/7.x/x86_64/pkgs/PASA/2.0.2
PASA=/rhome/cjinfeng/.pasa/

module load PASA/2.0.2
$PASA/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -R -g genome_sample.fasta -C -e 2 \
 -t all_transcripts.fasta.clean -T -u all_transcripts.fasta -f FL_accs.txt --ALIGNERS blat,gmap --CPU $cpu

module unload perl/5.20.2
module load perl/5.16.3
module load gmap/2014.12.28
module load kent/318
$PASA/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -R -g genome_sample.fasta -s 2 \
 -t all_transcripts.fasta.clean -T -u all_transcripts.fasta -f FL_accs.txt --ALIGNERS blat,gmap --CPU $cpu

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
