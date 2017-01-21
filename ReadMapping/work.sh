echo "map first 8 SMARTcell corrected reads"
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Csinensis_HZAU/Citrus_HZAU.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus.fastq --tool blasr --cpu 30 --output Citrus_sc8corrected_blasr --project Citrus_blasr --verbose > log 2>&1 &

echo "map first 24 SMARTcell corrected reads, subset"
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Csinensis_HZAU/Citrus_HZAU.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus.24sm_corrected_subset.fasta --tool blasr --cpu 30 --output Citrus_sc8corrected_blasr --project Citrus_24sm_blasr --verbose > log 2>&1 &

echo "yeast"
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_norm.fasta --tool blasr --cpu 30 --output yeast_norm_blasr --project yeast_norm -verbose > log 2>&1 &
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_4falcon.fasta --tool blasr --cpu 30 --output yeast_4falcon_blasr --project yeast_4falcon -verbose > log 2>&1 &
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_subread.fasta --tool blasr --cpu 30 --output yeast_subread_blasr --project yeast_subread -verbose > log 2>&1 &
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_MHAPcor.fastq --tool blasr --cpu 30 --output yeast_MHAP_blasr --project yeast_MHAP -verbose > log 2>&1 &
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_ass.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_4falcon.fasta --tool blasr --cpu 30 --project yeast_4falcon_ass -verbose > log 2>&1 &
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_ass.fa --tool blasr --cpu 30 --split 100 --project yeast_ass2w303 --genomealign -verbose > log 2>&1 &
echo "yeast quiver"
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_ass_quiver_round1.fasta --tool blasr --cpu 30 --split 100 --project yeast_ass_quiver_round1_2w303 --genomealign -verbose > log 2>&1 &

echo "Citrus_sm40"
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_asm40_quiver_round1.fasta --tool blasr --cpu 60 --split 100 --genomealign --project Citrus_asm40_quiver_round1 --verbose
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_asm40_quiver_round2.fasta --tool blasr --cpu 60 --split 100 --genomealign --project Citrus_asm40_quiver_round2 --verbose
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/citrus_ass_v1_line50_clean.fa --tool blasr --cpu 60 --split 100 --genomealign --project Citrus_asm40_raw --verbose
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/preads4falcon.fasta --tool blasr --cpu 60 --project Citrus_asm40_preads --verbose

#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_hybrid_scaffold.fa --tool blasr --cpu 30 --project Citrus_hybrid_scaffold --verbose
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_runCA.fa --split 700 --tool blasr --cpu 30 --project Citrus_runCA --verbose
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_hybrid_scaffold.fa --split 100 --tool blasr --cpu 30 --project Citrus_hybrid_scaffold --verbose
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_canu.fa --split 100 --tool blasr --cpu 30 --project Citrus_canu --verbose

perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_asm40_quiver_round2_merge.fasta --tool bwa --cpu 30 --split 2 --project Citrus_asm40_quiver_round2_merge --verbose
#merge but not remove tandem duplicate or break misjoin 
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_asm40_quiver_round2_merge_only.fasta --tool bwa --cpu 30 --split 2 --project Citrus_asm40_quiver_round2_merge_only --verbose > log 2>&1 &
#merge break level =1 Citrus_asm40_quiver_round2_merge_only_nobreak.fasta
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_asm40_quiver_round2_merge_only_nobreak.fasta --tool bwa --cpu 30 --split 2 --project Citrus_asm40_quiver_round2_merge_only_nobreak --verbose > log 2>&1 &

echo "PBcR assembly"
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/citrus_PBcR_v1_line50.fa --tool bwa --cpu 30 --split 2 --project citrus_PBcR_v1 --verbose > log 2>&1 &


echo "coverage along chromosome"
samtools view -bq 2 Citrus_asm40_quiver_round2.bam > Citrus_asm40_quiver_round2.filter_MQ2.bam &
bedtools makewindows -g ../Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.9.chrlen -w 1000000 -s 200000 > Cclementina_v1.0_scaffolds.9.slidingwin.bed
bedtools coverage -abam Citrus_asm40_quiver_round2.filter_MQ2.bam -b Cclementina_v1.0_scaffolds.9.slidingwin.bed -hist > Cclementina_v1.0_scaffolds.9.slidingwin.bed.hist.cov &
python Bam_Cov_SlidingWin.py --input Cclementina_v1.0_scaffolds.9.slidingwin.bed.hist.cov
python Convert_position_to_one_chromosome_bed.py --chr ../Citrus_diversity/bin/Cclementina_v1.0_scaffolds.chrlen --bed Cclementina_v1.0_scaffolds.9.slidingwin.bed.hist.cov.bed > Cclementina_v1.0_scaffolds.9.slidingwin.bed.hist.cov.1chr.bed
cat Cclementina_v1.0_scaffolds.9.slidingwin.bed.hist.cov.1chr.R | R --slave

echo "Neurospora"
#perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Neurospora_10/FungiDB-26_Ncrassa_OR74A_Genome.fasta -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON/Neurospora_ass/OR_falcon/falcon_mix3_1000_1000/p_ctg.fa --tool blasr --cpu 30 --split 100 --genomealign --project Nc_mix3 --verbose
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Neurospora_10/FungiDB-26_Ncrassa_OR74A_Genome.fasta -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON/Neurospora_ass/OR_falcon/data/R079-F02.1-Reads.line50.fasta --tool blasr --cpu 30 --split 2000 --genomealign --project Nc_OR_reads

echo "human, na12878" na12878_hybrid_scaffold.fa
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/hg18/hg18.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/na12878_hybrid_scaffold.fa --tool blasr --cpu 30 --split 10 --genomealign --project na12878_hybrid_scaffold --verbose
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/hg18/hg18.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/na12878_raw.fa --tool blasr --cpu 30 --split 40 --genomealign --project na12878_raw --verbose
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/hg18/hg18.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/na12878_falcon.fa --tool blasr --cpu 30 --split 40 --genomealign --project na12878_falcon --verbose

echo "citrus_ass_v3"
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/citrus_ass_v3_line50_clean.fa --tool bwa --cpu 30 --split 2 --project citrus_ass_v3 --verbose > log 2>&1 &

echo "citrus_canu1.3"
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/citrus.contigs.fasta --tool bwa --cpu 30 --split 2 --project citrus_canu1_3_ass --verbose > log 2>&1 &

echo "preads bwa mem for SNPsplit"
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/preads4falcon.fasta --tool bwa --cpu 30 --split 200 --project citrus_preads_bwa_mem --verbose > log 2>&1 &


echo "Done"
