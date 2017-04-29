echo "scaffold_2 testing"
samtools view ~/BigData/00.RD/Assembly/Pacbio/ReadMapping/Citrus_asm40_quiver_round2.bam | awk '$3~/scaffold_2$/ && $5 > 10' | cut -f1 | sed 's/\/.*//' > scaf2
perl ~/BigData/software/bin/getidseq.pl -l scaf2 -f citrus_ass_v1_quiver_round2.fasta.RepeatMasker.masked -o citrus_ass_v1_quiver_round2.fasta.RepeatMasker.scaf2.masked

echo "20kb"
ln -s Fairchild_falconv3_20kb_cov2_p_ctg_quiver_round1_pilon.fasta.RepeatMasker.masked.gz genome.fa.gz

