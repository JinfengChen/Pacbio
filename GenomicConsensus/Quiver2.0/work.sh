python Run_Quiver.py --bamdir merge_bam_split_contig --ref citrus_PBcR_v1_line50_split_contig --project Citrus_PBcR_quiver > log 2>&1 &
perl Finish_Quiver.pl --quiver Citrus_PBcR_quiver.consensus.fasta --falcon citrus_PBcR_v1.fa --output citrus_PBcR_v1_quiver_round1.fasta

