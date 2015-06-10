#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l mem=140gb
#PBS -l walltime=200:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

module load java/1.8.0_25

PATH=$PATH:/bigdata/cjinfeng/00.RD/Assembly/Pacbio/install/FALCON/DALIGNER/DALIGNER:/bigdata/cjinfeng/00.RD/Assembly/Pacbio/install/FALCON/DAZZ_DB/DAZZ_DB

python fc_run.py fc_run_ecoli.cfg 

#ref=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa
#asm=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON/yeast_test/2-asm-falcon/p_ctg.fa
#/usr/local/bin/dnadiff $ref $asm

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

