#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=50:00:00
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
#cpu=12 # need to between 1 and 16
#module load PASA/2.0.2
#seqclean transcripts.fasta -c $cpu
#PASA=/opt/linux/centos/7.x/x86_64/pkgs/PASA/2.0.2
PASA=/rhome/cjinfeng/.pasa/
#genome=Fairchild_v1.fasta
genome=Fairchild_v1.chr4.fasta
transcript=transcripts.fasta

module load PASA/2.0.2
$PASA/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -R -g $genome -C -e 2 \
 -t $transcript\.clean -T -u $transcript --ALIGNERS blat,gmap --CPU $CPU

module unload perl/5.20.2
module load perl/5.16.3
module load gmap/2014.12.28
module load kent/318
$PASA/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -R -g $genome -s 2 \
 -t $transcript\.clean -T -u $transcript --ALIGNERS blat,gmap --CPU $CPU

pasadb=cjinfeng
$PASA/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta ${pasadb}.assemblies.fasta --pasa_transcripts_gff3 ${pasadb}.pasa_assemblies.gff3

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
