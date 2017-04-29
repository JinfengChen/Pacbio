#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=60G
#SBATCH --time=40:00:00
#SBATCH -p highmem
#SBATCH --workdir=./
#SBATCH --output=splitSNP.stdout

#cd $PBS_O_WORKDIR
#qsub -t 1-7 indexbam_qsub.sh
start=`date +%s`

#blasr --nproc $PBS_NP test.fasta  --bestn 1 -m 5 --minMatch 19 --out
#python `pwd`/splitSNP/splitSNP_pacbio.py citrus_canu1_3_ass.bam t1 scaffold_1:47220:C:T

export PYTHONPATH=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/pythonlib/:/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/pythonlib/lib64/python2.7/site-packages/:/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/pythonlib/lib/python2.7/site-packages/:$PYTHONPATH

#python splitSNP/splitSNP_pipe.py --input scaffold_all.snp.list --bam citrus_20kb_bwa_mem.bam --size 5000 --output Pacbio_20kb_haplotype_reads
python splitSNP/splitSNP_extract_reads_Pacbio.py --input Pacbio_20kb_haplotype_reads --fasta citrus.pacbio_20kb.fasta
#python splitSNP/splitSNP_extract_reads_Pacbio.py --input Pacbio_20kb_haplotype_reads --fasta_list input.20kb.fofn

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

