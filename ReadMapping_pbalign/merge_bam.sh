#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=40gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

export LD_LIBRARY_PATH=/bigdata/stajichlab/cjinfeng/00.RD/Assembly/Pacbio/install/pitchfork/deployment/lib:$LD_LIBRARY_PATH
export PATH=/bigdata/stajichlab/cjinfeng/00.RD/Assembly/Pacbio/install/pitchfork/deployment/bin:$PATH

ls `pwd`/Citrus_PBcR_round1/Citrus_PBcR_round1.m150*_p0.bam > Citrus_PBcR_round1.bam.fofn
pbmerge -o Citrus_PBcR_round1.merged.bam Citrus_PBcR_round1.bam.fofn
samtools index Citrus_PBcR_round1.merged.bam
samtools view -H Citrus_PBcR_round1.merged.bam | awk '{print $2}' | grep "^SN" | awk '{gsub(/SN\:/,""); print}' > Citrus_PBcR_round1.merged.contigs.txt
if [ ! -d Citrus_PBcR_round1.merged.bam_split_contig ]; then
    mkdir Citrus_PBcR_round1.merged.bam_split_contig
fi

if [ ! -d citrus_PBcR_v1_line50_split_contig ]; then
    mkdir citrus_PBcR_v1_line50_split_contig
fi


for c in `cat Citrus_PBcR_round1.merged.contigs.txt` ; do
    echo processing $c
    samtools view -@ $PBS_NP -bh Citrus_PBcR_round1.merged.bam $c > Citrus_PBcR_round1.merged.bam_split_contig/$c.bam
    pbindex Citrus_PBcR_round1.merged.bam_split_contig/$c.bam
    samtools index Citrus_PBcR_round1.merged.bam_split_contig$c.bam
    perl ~/BigData/software/bin/fastaDeal.pl --get_id $c citrus_PBcR_v1_line50.fa > citrus_PBcR_v1_line50_split_contig/$c\.fa
    samtools faidx citrus_PBcR_v1_line50_split_contig/$c\.fa 
done


#dataset create --type AlignmentSet merge_bam.XML merge_bam.fofn
#dataset split --contigs --outdir merge_bam_split_by_contig merge_bam.XML



end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

