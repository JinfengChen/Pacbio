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

