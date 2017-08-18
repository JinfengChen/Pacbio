#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=add_locus_tag.sh.stdout
#SBATCH --job-name="add_locus_tag"
#SBATCH -p intel

module load bedops/2.4.24

genome=Fairchild_v1.fasta
noTE_pre=Fairchild.optimized_model.noTE_highqual_AS_best
locus_tag=BTW09_

awk '$3=="gene"' $noTE_pre\.gff > $noTE_pre\.gene.gff
gff2bed < $noTE_pre\.gene.gff > $noTE_pre\.gene.bed
python -m jcvi.annotation.reformat rename --prefix=$locus_tag $noTE_pre\.gene.bed 
mv BTW09_.ids $noTE_pre\.gene.bed.ids
python update_locus_id.py --gff $noTE_pre\.gff --id $noTE_pre\.gene.bed.ids --genome $genome

module load augustus
lineage=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/dataset/embryophyta_odb9
export AUGUSTUS_CONFIG_PATH=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/augustus_config
temp=/rhome/cjinfeng/Rice/Rice_population_sequence/BUSCO/PlantGenome/tmp
#protein=$OUTPUT\.noTE.pep.fa
#protein=$OUTPUT\.noTE_highqual_AS_best.pep.fa
protein=$noTE_pre\.locus_tag_id.pep.fa
output=$protein\.BUSCO

if [ -e $protein ] && [ ! -e "run_$output" ]; then
   python ~/BigData/00.RD/Assembly/Pacbio/install/BUSCO/busco/scripts/run_BUSCO.py --in $protein --cpu $SLURM_NTASKS --out $output --lineage_path $lineage --mode prot --tmp_path $temp
fi

echo "Done"


echo "Done"
