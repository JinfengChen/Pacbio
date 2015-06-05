#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR


start=`date +%s`

module load java/1.8.0_25

PATH=$PATH:/bigdata/cjinfeng/00.RD/Assembly/Pacbio/install/FALCON/DALIGNER/DALIGNER:/bigdata/cjinfeng/00.RD/Assembly/Pacbio/install/FALCON/DAZZ_DB/DAZZ_DB

python fc_run.py fc_run_ecoli.cfg 

#ref=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Ecoli/selfSampleData/reference.fasta
#asm=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/PBcR/K12/9-terminator/asm.ctg.fasta
#/usr/local/bin/dnadiff $ref $asm

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

