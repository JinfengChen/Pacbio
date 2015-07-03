#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR

start=`date +%s`

#/usr/local/bin/nucmer -mumreference -c 1000 -l 100 -d 10 W303.ref.fa yeast_ass.fa
#/usr/local/bin/delta-filter -1 out.delta > out.1delta
#/usr/local/bin/mummerplot -large -fat --postscript out.1delta
#set size 3,3; size not right in out.gp; remove or modify
#gnuplot out.pg

module load mummer
mummer -mum -l 200 -b -c Fairchild.fasta Cclementina_v1.0_scaffolds.fa > Fairchild_JGI.mummer.mums

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

