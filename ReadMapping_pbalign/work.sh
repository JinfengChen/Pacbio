echo "yeast"
ls /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/Pacbio_assembly/*/Analysis_Results/*.bax.h5 > input.yeast.fofn
perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/yeast_ass.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.yeast.fofn --output yeast_ass_reads --project yeast_ass_reads > log 2>&1 &


echo "Citrus"
perl step1_Mapping_h5_pbalign.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/Citrus_24sm_ass.fa --input /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping_pbalign/input.fofn --output Citrus_24sm_ass --project Citrus_24sm_reads > log 2>&1 &

