#align
cat ~/BigData/00.RD/Assembly/Pacbio/GenomeAlign_mummer/assembly_falconv3_20kb_cov2_p_ctg_quiver_round1_pilon_haplomerge_10xbionabo_hybrid/*.1align | grep -v "bigdata" | grep "chr" > citrus_FCM_CLM.1align
#gap
perl ~/BigData/software/bin/fasta2agp.pl ~/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.chr.fa > CLM.scaf.agp
sed 's/scaffold_/chr/' CLM.scaf.agp > CLM.chr.agp
perl ~/BigData/software/bin/fasta2agp.pl ~/BigData/05.Share/Fairchild/Fairchild_v1.0.chr1-9.fasta > FCM.chr.agp
#chrlen
head -n 9 ~/BigData/00.RD/Assembly/Pacbio/Reference/Fairchild/Fairchild_v1.fasta.len > FCM.chrlen
head -n 9 ~/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.chrlen > CLM.chrlen
sed 's/scaffold_/chr/' CLM.chrlen > CLM.chr.chrlen
perl drawchr.pl --align citrus_FCM_CLM.1align > log 2>&1 &

