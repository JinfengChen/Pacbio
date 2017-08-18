wget http://www.hrt.msu.edu/uploads/535/78637/Tpases020812.gz
wget http://www.hrt.msu.edu/uploads/535/78637/Tpases020812DNA.gz
gunzip Tpases020812.gz
gunzip Tpases020812DNA.gz
#HWB, gene=30122, mRNA=42886
#CSI, gene=29655, mRNA=44275
#CLM, gene=24533, mRNA=33929

echo "gene rename" 
module load bedops/2.4.24
awk '$3=="gene"' /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Annotation/gene/EVM/run_from_maker/run_FCM/weight_PASA_transdecoder_maker_genewise_OP/EVM.all.gff3 > FCM_all.pep.gene.gff
gff2bed < FCM_all.pep.gene.gff > FCM_all.pep.gene.bed
python -m jcvi.annotation.reformat rename --prefix=BTW09_ FCM_all.pep.gene.bed


gff2bed < test.gff3 > test.bed
python -m jcvi.annotation.reformat rename test.bed

echo "check completeness as AS analysis"
perl getGene.pl FCM_all.pep.noTE.AS_10kb.gff Fairchild_v1.fasta > FCM_all.pep.noTE.AS_10kb.cds.fa
perl ~/BigData/software/bin/cds2aa.pl FCM_all.pep.noTE.AS_10kb.cds.fa > FCM_all.pep.noTE.AS_10kb.pep.fa
ls FCM_all.pep.noTE.AS_10kb.pep.fa > genome.list
sbatch --array 1 buscov3_gene_array.sh

echo "best gene model"
sbatch filter_low_quality_gene.sh
#number of gene in clementine has AS: 4984 
awk '$3=="mRNA"' ~/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_gene.chr_edited.gff3 | grep "longest=0" | sed 's/.*Parent=//'| uniq | sort | uniq | wc -l
#numebr of gene in fairchild has AS: 7798 
#awk '$3=="mRNA"' Fairchild.optimized_model.noTE_highqual_AS.gff | cut -f9 | sed 's/scaffold//g' | grep "f" | sed 's/.*Parent=//' | uniq | sort | uniq | wc -l
#number of gene in fairchild has AS_best: 5953
awk '$3=="mRNA"' Fairchild.optimized_model.noTE_highqual_AS_best.gff | cut -f9 | sed 's/scaffold//g' | grep "f" | sed 's/.*Parent=//' | uniq | sort | uniq | wc -l

#manual selection
grep "^fusion: 4" filter_low_quality_gene.sh.stdout.gene_fusion -A 10 | less -S
grep "^fusion: 3" filter_low_quality_gene.sh.stdout.gene_fusion -A 10 | less -S
#check in jbrowse based on position
grep "replaced" filter_low_quality_gene.sh.stdout -B 2

#large intron
perl ~/BigData/software/bin/gff2intron.pl --gff Fairchild.optimized_model.noTE_highqual_AS_best.gff > Fairchild.optimized_model.noTE_highqual_AS_best.intron.gff
awk '$3=="intron" && $5-$4 > 10000' Fairchild.optimized_model.noTE_highqual_AS_best.intron.gff | less -S
awk '$3=="intron" && $5-$4 > 10000' Fairchild.optimized_model.noTE_highqual_AS_best.intron.gff | awk '{print $5-$4}' | less -S
python long_intron_gene.py --gff Fairchild.optimized_model.noTE_highqual_AS_best.gff > Fairchild.optimized_model.noTE_highqual_AS_best.gff.long_intron.gene.list &


#rename gff
sbatch add_locus_tag.sh

