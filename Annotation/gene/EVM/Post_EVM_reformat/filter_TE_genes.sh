#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00
#SBATCH --output=filter_TE_genes.sh.stdout
#SBATCH --job-name="filter_TE_genes"
#SBATCH -p intel

module load ncbi-blast/2.2.26

SWISS_PROT_PATH="Tpases020812"
if [ ! -e $SWISS_PROT_PATH\.phr ]; then    
   formatdb -i $SWISS_PROT_PATH
fi

pep=FCM_all.pep.fa
pre=${pep%.fa}
noTE=$pre\.noTE.fa
if [ ! -e $pep\_blast_results.txt ]; then
    blastall -p blastp -i $pep -d $SWISS_PROT_PATH -e 1e-10 -o $pep\_blast_results.txt -a $SLURM_NTASKS
fi

#perl ~/BigData/00.RD/Annotation/HEG4/protein-map-genome/bin/solar/solar.pl -d -1 $pep\_blast_results.txt > $pep\_blast_results.solar
#perl ~/BigData/software/bin/bestAlign.pl $pep\_blast_results.solar -cutoff 0.3 | cut -f1 > $pep\_blast_results.solar.id
#perl ~/BigData/software/bin/getidseq.pl -l $pep\_blast_results.solar.id -f $pep -o $noTE -r
perl ~/BigData/software/bin/fastaDeal.pl --attr id $noTE > $noTE\.id
python -m jcvi.formats.gff parents $pre\.gff $noTE\.id > $noTE\.mrna_gene.table
cut -f2 $noTE\.mrna_gene.table > $noTE\.gene.id
perl getidgff.pl -l $noTE\.gene.id -g $pre\.gff -o $pre\.noTE.gff

module load augustus
lineage=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/dataset/embryophyta_odb9
export AUGUSTUS_CONFIG_PATH=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/augustus_config
temp=/rhome/cjinfeng/Rice/Rice_population_sequence/BUSCO/PlantGenome/tmp
output=$pre\.BUSCO

if [ ! -e "run_$output" ]; then
   python ~/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/scripts/run_BUSCO.py --in $noTE --cpu $SLURM_NTASKS --out $output --lineage_path $lineage --mode prot --tmp_path $temp
fi
