#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=10:00:00
#SBATCH --output=add_locus_tag.sh.stdout
#SBATCH --job-name="add_locus_tag"
#SBATCH -p intel

module load bedops/2.4.24

noTE_pre=FCM_all.pep.noTE
locus_tag=BTW09_

awk '$3=="gene"' $noTE_pre\.gff > $noTE_pre\.gene.gff
gff2bed < $noTE_pre\.gene.gff > $noTE_pre\.gene.bed
python -m jcvi.annotation.reformat rename --prefix=$locus_tag $noTE_pre\.gene.bed 


echo "Done"
