#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

module load freebayes
freebayes -f ../../../Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -@ all.snp.vcf.gz WMM.bam > WMM.genotype.vcf

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

