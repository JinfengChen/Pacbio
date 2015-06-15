#!/bin/bash
#PBS -l nodes=1:ppn=32
#PBS -l mem=160gb
#PBS -l walltime=200:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

module load java/1.8.0_25

PBcR_bin=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/wgs-8.3rc2/Linux-amd64/bin
spec_file=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/PBcR/Citrus_4falcon.spec
lib_name=Citrus_4falcon
gsize=370000000

#fast consensus (-pbCNS parameter)
#$PBcR_bin/PBcR -threads $PBS_NP -length 500 -partitions 200 -l $lib_name -s $spec_file -fastq $fq_file genomeSize=$gsize
#default
#$PBcR_bin/runCA -s $spec_file -p asm -d $lib_name ovlRefBlockLength=100000000000 ovlRefBlockSize=0 useGrid=0 scriptOnGrid=0 unitigger=bogart ovlErrorRate=0.03 utgErrorRate=0.025 cgwErrorRate=0.1 cnsErrorRate=0.1 utgGraphErrorLimit=0 utgGraphErrorRate=0.025 utgMergeErrorLimit=0 utgMergeErrorRate=0.025 frgCorrBatchSize=100000 doOverlapBasedTrimming=1 obtErrorRate=0.03 obtErrorLimit=4.5 frgMinLen=3000 ovlMinLen=100 "batOptions=-RS -NS -CS" consensus=pbutgcns merSize=22 cnsMaxCoverage=1 cnsReuseUnitigs=1 gridEnginePropagateHold="pBcR_asm"  $lib_name.frg 
#merge haplotype
$PBcR_bin/runCA -s $spec_file -p asm -d $lib_name ovlRefBlockLength=100000000000 ovlRefBlockSize=0 useGrid=0 scriptOnGrid=0 unitigger=bogart ovlErrorRate=0.03 utgErrorRate=0.025 cgwErrorRate=0.1 cnsErrorRate=0.1 utgGraphErrorLimit=3.25 utgGraphErrorRate=0.07 utgMergeErrorLimit=5.25 utgMergeErrorRate=0.085 frgCorrBatchSize=100000 doOverlapBasedTrimming=0 obtErrorRate=0.03 obtErrorLimit=4.5 frgMinLen=3000 ovlMinLen=100 "batOptions=-RS -NS -CS" consensus=pbutgcns merSize=22 cnsMaxCoverage=1 cnsReuseUnitigs=1 gridEnginePropagateHold="pBcR_asm"  $lib_name.frg

#ref=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Ecoli/selfSampleData/reference.fasta
#asm=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/PBcR/K12/9-terminator/asm.ctg.fasta
#/usr/local/bin/dnadiff $ref $asm

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

