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

chr=7
prefix=Citrus_chr$chr
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/nucmer -mumreference --forward -c 1000 -l 100 -d 10 -p $prefix Cclementina_v1.0_scaffolds.chr$chr\.fasta Fairchild.chr$chr\.fasta
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/delta-filter -1 $prefix\.delta > $prefix\.1delta
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -fat --postscript -p $prefix\.clean $prefix\.1delta
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/mummerplot -fat --postscript -p $prefix $prefix\.delta

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

