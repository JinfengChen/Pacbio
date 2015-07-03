echo "test and debug quiver"
#1. ref fasta need to 50 or some other base per line, not single long line per sequence
#2. ref fasta name need to short and do not have space
#3. do not use one SMRTcell to merge and sort
#4. --mapQvThreshold=0 in denovo assembly consensus to deal with low quality read mapping in repeat regions
quiver -j1 aligned_reads.cmp.h5 -p C2.NoQVsModel -r HCV_Ref_For_187140.fasta -o test.gff -o test.fasta -o test.fastq
python test_ref.py --input yeast_ass_reads.cmp.h5 --ref yeast_ass_line50.fa > log 2>&1 &
python Run_Quiver.py --cmpdir ../../ReadMapping_pbalign/yeast_ass_reads_cmp_split --ref yeast_ass_line50_clean.fa --project yeast_quiver > yeast_quiver.log 2>&1 &

echo "check info in cmp.h5"
h5ls -r 000028F.sort.cmp.h5 | grep MovieInfo


echo "Citrus_ass36_test3"
ln -s ~/BigData/00.RD/Assembly/Pacbio/FALCON/Citrus_test/2-asm-falcon_test3/p_ctg.fa citrus_ass.fa
perl ~/software/bin/fastaDeal.pl --attr id citrus_ass.fa > citrus_ass.id
perl ~/software/bin/fastaDeal.pl --reform line50 citrus_ass.fa > citrus_ass_line50.fa
perl getidseq.pl -l citrus_ass.id -f citrus_ass_line50.fa -o citrus_ass_line50_clean.fa
samtools faidx citrus_ass_line50_clean.fa

echo "Citrus_ass40_test3"
ln -s /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON/Citrus_V1/Citrus_v1_raw.fasta citrus_ass_v1.fa
perl ~/software/bin/fastaDeal.pl --attr id citrus_ass_v1.fa > citrus_ass_v1.id
perl ~/software/bin/fastaDeal.pl --reform line50 citrus_ass_v1.fa > citrus_ass_v1_line50.fa
perl getidseq.pl -l citrus_ass.id -f citrus_ass_v1_line50.fa -o citrus_ass_v1_line50_clean.fa
samtools faidx citrus_ass_v1_line50_clean.fa
python Run_Quiver.py --cmpdir ../../ReadMapping_pbalign/Citrus_40sm_reads_cmp_split/ --ref citrus_ass_v1_line50_clean.fa --project Citrus_asm40_quiver > log 2>&1 &

