echo "yeast"
#ls /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/Pacbio_assembly/*/Analysis_Results/*.bax.h5 > input.yeast.fofn
#perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/yeast_ass_line50_clean.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.yeast.fofn --output yeast_ass_reads --project yeast_ass_reads --step 3 --verbose > log 2>&1 &

perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/yeast_ass_line50_clean.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.yeast_1.fofn --output yeast_ass_reads --project yeast_ass_reads --step 12 --verbose > log 2>&1 &

echo "Citrus_asm36"
#ls /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Citrus/SMARTcell_raw/m*/Analysis_Results/*.bax.h5 > input.Citrus.fofn
#perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/citrus_ass_line50_clean.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.Citrus.fofn --project Citrus_36sm_reads --step 3 --verbose  > log 2>&1 &

echo "Citrus_asm40"
#round1
#perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/citrus_ass_v1_line50_clean.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.Citrus_highqual.fofn --project Citrus_40sm_reads --step 23 --verbose  > log 2>&1 &
#round2
#perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/citrus_ass_v1_quiver_round1.fasta --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.Citrus_highqual.fofn --project Citrus_40sm_reads_round2_1 --step 12 --verbose  > log 2>&1 &


echo "Citrus_PBcR"
perl step1_Mapping_h5_pbalign_bam.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/citrus_PBcR_v1_line50.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.Citrus_A02_1.fofn --project Citrus_PBcR_round1 --step 12 --verbose  > log 2>&1 
#
perl step1_Mapping_h5_pbalign_bam.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/citrus_PBcR_v1_line50.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.Citrus_highqual.fofn --project Citrus_PBcR_round1 --step 12 --verbose  > log 2>&1
