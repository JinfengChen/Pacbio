echo "map first 8 SMARTcell corrected reads"
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Csinensis_HZAU/Citrus_HZAU.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus.fastq --tools blasr --cpu 30 --output Citrus_sc8corrected_blasr --project Citrus_blasr --verbose > log 2>&1 &

echo "map first 24 SMARTcell corrected reads, subset"
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Csinensis_HZAU/Citrus_HZAU.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus.24sm_corrected_subset.fasta --tools blasr --cpu 30 --output Citrus_sc8corrected_blasr --project Citrus_24sm_blasr --verbose > log 2>&1 &

echo "yeast"
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_norm.fasta --tools blasr --cpu 30 --output yeast_norm_blasr --project yeast_norm -verbose > log 2>&1 &
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_4falcon.fasta --tools blasr --cpu 30 --output yeast_4falcon_blasr --project yeast_4falcon -verbose > log 2>&1 &
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_subread.fasta --tools blasr --cpu 30 --output yeast_subread_blasr --project yeast_subread -verbose > log 2>&1 &
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_MHAPcor.fastq --tools blasr --cpu 30 --output yeast_MHAP_blasr --project yeast_MHAP -verbose > log 2>&1 &
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_ass.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_4falcon.fasta --tools blasr --cpu 30 --project yeast_4falcon_ass -verbose > log 2>&1 &
perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Yeast/W303.ref.fa -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/ReadMapping/yeast_ass.fa --tools blasr --cpu 30 --split 100 --project yeast_ass2w303 --genomealign -verbose > log 2>&1 &

