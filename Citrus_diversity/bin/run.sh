#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

strain=PKM
module load vcftools
zcat ../input/SNP/$strain\.vcf.gz| vcf-to-tab > $strain\.vcf.tab
python HeterozygousDensity.py --input $strain\.vcf.tab > $strain\.vcf.tab.het.density
python Convert_position_to_one_chromosome_tab.py --chr Cclementina_v1.0_scaffolds.chrlen --tab $strain\.vcf.tab.het.density > $strain\.vcf.tab.het.1chr.density
python plot_chr_density.py --input $strain\.vcf.tab.het.1chr.density

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

