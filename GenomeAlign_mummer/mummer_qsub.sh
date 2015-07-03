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

prefix=Citrus_sm40
#nucmer: slow but align with reference order: takes 40min for 370 Mb genome, that's okay.
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/nucmer -mumreference -c 1000 -l 200 -d 10 -p $prefix Cclementina_v1.0_scaffolds.fa Fairchild.fasta
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/delta-filter -1 $prefix\.delta > $prefix\.1delta
#mummer: fast for large genome, but order is ranked by name and too much signal, need to filter. okay for chromosome level
#/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummer -mumreference -l 200 -b -c Cclementina_v1.0_scaffolds.fa Fairchild.fasta > $prefix\.mummer


#/usr/local/bin/delta-filter -1 out.delta > out.1delta
#/usr/local/bin/mummerplot -large -fat --postscript out.1delta
#set size 1,1; size not right in out.gp; remove or modify
#gnuplot out.pg
#nucmer works with contig and reference comparsion, which is multi to multi
#mummer works with 1 to 1 comparision

#raw plot
#/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -R Cclementina_v1.0_scaffolds.chr -Q Fairchild.chr -fat --postscript -p $prefix $prefix\.delta
#/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -R Cclementina_v1.0_scaffolds.chr -Q Fairchild.scaf -fat --postscript -p $prefix\.1 $prefix\.delta
#/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -R Cclementina_v1.0_scaffolds.scaf -Q Fairchild.chr -fat --postscript -p $prefix\.2 $prefix\.delta
#clean plot: -fat will make strand confuse, -large will make size problem
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -R Cclementina_v1.0_scaffolds.chr -Q Fairchild.chr --postscript -p $prefix\.clean $prefix\.1delta
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -R Cclementina_v1.0_scaffolds.chr -Q Fairchild.scaf --postscript -p $prefix\.clean.1 $prefix\.1delta
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -R Cclementina_v1.0_scaffolds.scaf -Q Fairchild.chr --postscript -p $prefix\.clean.2 $prefix\.1delta
#mumer plot
#/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -R Cclementina_v1.0_scaffolds.fa -Q Fairchild.fasta --layout -fat --postscript -p $prefix\.mummer $prefix\.mummer

#chr vs chr plot
for i in {1..9}; 
do
    /rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -r scaffold_$i -q chr$i --postscript -p $prefix\.chr$i $prefix\.delta
    /rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -r scaffold_$i -q chr$i --postscript -p $prefix\.clean.chr$i $prefix\.1delta
done

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

