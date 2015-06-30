echo "yeast"
#ls /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/Pacbio_assembly/*/Analysis_Results/*.bax.h5 > input.yeast.fofn
#perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/yeast_ass_line50_clean.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.yeast.fofn --output yeast_ass_reads --project yeast_ass_reads --step 3 --verbose > log 2>&1 &

echo "Citrus"
#ls /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Citrus/SMARTcell_raw/m*/Analysis_Results/*.bax.h5 > input.Citrus.fofn
perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/citrus_ass_line50_clean.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.Citrus.fofn --project Citrus_36sm_reads --step 123 --verbose  > log 2>&1 &

