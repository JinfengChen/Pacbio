#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=4G
#SBATCH --time=40:00:00
#SBATCH --output=filter_low_quality_gene.sh.stdout
#SBATCH --job-name="filter_low_quality_gene"
#SBATCH -p intel

module load ncbi-blast/2.2.26

SWISS_PROT_PATH="uniprot_sprot.fasta"
if [ ! -e $SWISS_PROT_PATH\.phr ]; then    
   formatdb -i $SWISS_PROT_PATH
fi

TPASE_PATH="Tpases020812"
if [ ! -e $TPASE_PATH\.phr ]; then    
   formatdb -i $TPASE_PATH
fi


pep=FCM_all.pep.fa
pre=${pep%.fa}
if [ ! -e $pep\_uniprot_sprot_blast_results.txt ]; then
    blastall -p blastp -i $pep -d $SWISS_PROT_PATH -e 1e-10 -o $pep\_uniprot_sprot_blast_results.txt -a $SLURM_NTASKS
fi

if [ ! -e $pep\_Tpase_blast_results.txt ]; then
    blastall -p blastp -i $pep -d $TPASE_PATH -e 1e-10 -o $pep\_Tpase_blast_results.txt -a $SLURM_NTASKS
fi

#EVM_GFF=FCM_all.pep.gff
GENOME=Fairchild_v1.fasta
EVM_GFF=test.gff
#MAKER_GFF=Fairchildv1.all.gff
MAKER_GFF=ABINITIO_PREDICTION.gff
PASA_AS_GFF=FCM_all.pep.noTE.AS_10kb.gff
REPEAT_GFF=Fairchild_v1.fasta.RepeatMasker.out.gff
TE_BLASTP=$pep\_Tpase_blast_results.txt
US_BLASTP=$pep\_uniprot_sprot_blast_results.txt
OUTPUT=Fairchild.optimized_model


python optimize_gene_model.py --evm_gff $EVM_GFF \
                              --maker_gff $MAKER_GFF \
                              --pasa_as_gff $PASA_AS_GFF \
                              --repeat_gff $REPEAT_GFF \
                              --TE_blastp $TE_BLASTP \
                              --US_blastp $US_BLASTP \
                              --genome $GENOME \
                              --output $OUTPUT 


#module unload perl
#module load maker/2.31.8
#cp $pre\.gff $pre\.rename.gff
#maker_functional_gff $SWISS_PROT_PATH $pep\_uniprot_sprot_blast_results.txt $pre\.rename.gff > $pre\.rename.putative_function.gff

#perl ~/BigData/00.RD/Annotation/HEG4/protein-map-genome/bin/solar/solar.pl -d -1 $pep\_blast_results.txt > $pep\_blast_results.solar
#perl ~/BigData/software/bin/bestAlign.pl $pep\_blast_results.solar -cutoff 0.3 | cut -f1 > $pep\_blast_results.solar.id
#perl ~/BigData/software/bin/getidseq.pl -l $pep\_blast_results.solar.id -f $pep -o $noTE -r
#perl ~/BigData/software/bin/fastaDeal.pl --attr id $noTE > $noTE\.id
#python -m jcvi.formats.gff parents $pre\.gff $noTE\.id > $noTE\.mrna_gene.table
#cut -f2 $noTE\.mrna_gene.table > $noTE\.gene.id
#perl getidgff.pl -l $noTE\.gene.id -g $pre\.gff -o $pre\.noTE.gff

module load augustus
lineage=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/dataset/embryophyta_odb9
export AUGUSTUS_CONFIG_PATH=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/augustus_config
temp=/rhome/cjinfeng/Rice/Rice_population_sequence/BUSCO/PlantGenome/tmp
#protein=$OUTPUT\.noTE.pep.fa
protein=$OUTPUT\.noTE_highqual_AS_best.pep.fa
output=$protein\.BUSCO

if [ -e $protein ] && [ ! -e "run_$output" ]; then
   python ~/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/scripts/run_BUSCO.py --in $protein --cpu $SLURM_NTASKS --out $output --lineage_path $lineage --mode prot --tmp_path $temp
fi

echo "Done"
