echo "Citrus_sm40"
#whole genome
perl ~/software/bin/fastaDeal.pl --attr id:len Fairchild.fasta | grep "chr" | awk '{print $1"\t"$2"\t+"}' > Fairchild.chr
perl ~/software/bin/fastaDeal.pl --attr id:len Fairchild.fasta | grep "chr" -v | awk '{print $1"\t"$2"\t+"}' > Fairchild.scaf
perl ~/software/bin/fastaDeal.pl --attr id:len Cclementina_v1.0_scaffolds.fa | head -n 9 | awk '{print $1"\t"$2"\t+"}' > Cclementina_v1.0_scaffolds.chr
perl ~/software/bin/fastaDeal.pl --attr id:len Cclementina_v1.0_scaffolds.fa | tail -n 1389 | awk '{print $1"\t"$2"\t+"}' > Cclementina_v1.0_scaffolds.scaf

qsub mummer_qsub.sh
#longest contig
perl ~/software/bin/fastaDeal.pl --pat 000000F ~/BigData/00.RD/Assembly/Pacbio/ALLMAPS/Citrus_asm40/scaffolds.fasta > 000000F.fasta
perl ~/software/bin/fastaDeal.pl --pat scaffold_3 Cclementina_v1.0_scaffolds.fa > chr3.fasta
qsub mummer_longest.sh
perl ~/software/bin/fastaDeal.pl --substruct 44188106-49301823 chr3.fasta > chr3_44188106_49301823.fasta
/rhome/cjinfeng/BigData/software/mummer/MUMmer3.23/dnadiff -p longest chr3_44188106_49301823.fasta 000000F.fasta

#single chromosome
perl ~/software/bin/fastaDeal.pl --get_id scaffold_2 Cclementina_v1.0_scaffolds.fa > Cclementina_v1.0_scaffolds.chr2.fasta
perl ~/software/bin/fastaDeal.pl --get_id scaffold_4 Cclementina_v1.0_scaffolds.fa > Cclementina_v1.0_scaffolds.chr4.fasta
perl ~/software/bin/fastaDeal.pl --get_id scaffold_7 Cclementina_v1.0_scaffolds.fa > Cclementina_v1.0_scaffolds.chr7.fasta
perl ~/software/bin/fastaDeal.pl --get_id chr2 Fairchild.fasta > Fairchild.chr2.fasta
perl ~/software/bin/fastaDeal.pl --get_id chr4 Fairchild.fasta > Fairchild.chr4.fasta
perl ~/software/bin/fastaDeal.pl --get_id chr7 Fairchild.fasta > Fairchild.chr7.fasta
perl ~/software/bin/fastaDeal.pl --reverse --complement Fairchild.chr2.fasta > Fairchild.chr2_rec.fasta
perl ~/software/bin/fastaDeal.pl --reverse --complement Fairchild.chr4.fasta > Fairchild.chr4_rec.fasta
perl ~/software/bin/fastaDeal.pl --reverse --complement Fairchild.chr7.fasta > Fairchild.chr7_rec.fasta


