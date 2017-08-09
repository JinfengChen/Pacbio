#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=10:00:00
#SBATCH --output=add_locus_tag.sh.stdout
#SBATCH --job-name="add_locus_tag"
#SBATCH -p intel

module load bedops/2.4.24

noTE_pre=Fairchild.optimized_model.noTE_highqual_AS_best
locus_tag=BTW09_

#awk '$3=="gene"' $noTE_pre\.gff > $noTE_pre\.gene.gff
#gff2bed < $noTE_pre\.gene.gff > $noTE_pre\.gene.bed
#python -m jcvi.annotation.reformat rename --prefix=$locus_tag $noTE_pre\.gene.bed 
#mv BTW09_.ids $noTE_pre\.gene.bed.ids
python update_locus_id.py --gff $noTE_pre\.gff --id $noTE_pre\.gene.bed.ids

echo "Done"
