#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR

start=`date +%s`

export PERL5LIB=/opt/linux/centos/7.x/x86_64/pkgs/mummer/3.23/scripts/:$PERL5LIB
module load mummer

#/usr/local/bin/nucmer -mumreference -c 1000 -l 100 -d 10 W303.ref.fa yeast_ass.fa
#/usr/local/bin/delta-filter -1 out.delta > out.1delta
#/usr/local/bin/mummerplot -large -fat --postscript out.1delta
#set size 3,3; size not right in out.gp; remove or modify
#gnuplot out.pg

mummer -mum -l 200 -b -c 000000F.fasta chr3.fasta > longest.mummer.mums
mummerplot -x "[0,5086785]" -y "[44188106,49301823]" --postscript -p longest longest.mummer.mums

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

