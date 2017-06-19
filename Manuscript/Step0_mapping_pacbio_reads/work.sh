perl step1_Mapping_large_blasr.pl --ref /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Reference/Fairchild/Fairchild_v1.fasta -1 /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Manuscript/Step0_mapping_pacbio_reads/preads4falcon.fasta --tool bwa --cpu 30 --split 5000 --project citrus_raw_plus_20kb_pread_bwa_mem --verbose > log 2>&1 &

