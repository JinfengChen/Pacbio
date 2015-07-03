#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=200gb
#PBS -l walltime=200:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR

#ppn=32;mem=240gb
#ppn=48;mem=400gb
#ppn=64;mem=500gb

start=`date +%s`

#module load java/8u25
#PATH=$PATH:/bigdata/stajichlab/cjinfeng/00.RD/Assembly/Pacbio/install/FALCON/DALIGNER/DALIGNER:/bigdata/stajichlab/cjinfeng/00.RD/Assembly/Pacbio/install/FALCON/DAZZ_DB/DAZZ_DB


cpu=$PBS_NP
#cpu=1
quiver=/opt/linux/centos/7.x/x86_64/pkgs/python/2.7.5/bin/quiver 
#quiver=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/GenomicConsensus/bin/quiver
#cmp=yeast_ass_reads.sort.cmp.h5
cmp=000028F.sort.cmp.h5
ref=yeast_ass_line50_clean.fa
prefix=yeast

#echo $cpu, $prefix
#$quiver -j$cpu $cmp -p C2.NoQVsModel -r $ref --referenceChunkSize 100000 --referenceChunkOverlap 20000 -o $prefix.gff -o $prefix.consensus.fasta -o $prefix.consensus.fastq
#$quiver -j$cpu $cmp -p C2.NoQVsModel -r $ref -o $prefix.gff -o $prefix.consensus.fasta -o $prefix.consensus.fastq
$quiver -j$cpu $cmp -r $ref -o $prefix.gff -o $prefix.consensus.fasta -o $prefix.consensus.fastq

#ref=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa
#asm=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON/yeast_test/2-asm-falcon/p_ctg.fa
#/usr/local/bin/dnadiff $ref $asm

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

