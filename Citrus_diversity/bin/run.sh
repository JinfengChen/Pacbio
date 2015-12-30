#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

module load vcftools
#files=(FAIRCHILD LAP CHP SSO SWO PKM CLM WLM WMM)
files=(CLM)
for strain in ${files[@]}; do
    #freebayes
    #zcat ../input/SNP/$strain\.vcf.gz| vcf-to-tab > $strain\.vcf.tab
    #python HeterozygousDensity.py --input $strain\.vcf.tab > $strain\.vcf.tab.het.density
    #python Convert_position_to_one_chromosome_tab.py --chr Cclementina_v1.0_scaffolds.chrlen --tab $strain\.vcf.tab.het.density > $strain\.vcf.tab.het.1chr.density
    #python plot_chr_density.py --input $strain\.vcf.tab.het.1chr.density
    #python AncesterAnalysis.py --input $strain\.vcf.tab > $strain\.vcf.tab.ancester.bed
    
    #GATK
    if [ ! -e Citrus_6strain.vcf.tab ]; then
       cat ../input/SNP/Citrus_6strain.vcf | vcf-to-tab > Citrus_6strain.vcf.tab
    fi
    python HeterozygousDensity.py --input Citrus_6strain.vcf.tab --strain $strain > $strain\.vcf.tab.het.GATK.density
    python Convert_position_to_one_chromosome_tab.py --chr Cclementina_v1.0_scaffolds.chrlen --tab $strain\.vcf.tab.het.GATK.density > $strain\.vcf.tab.het.1chr.GATK.density
    python plot_chr_density.py --input $strain\.vcf.tab.het.1chr.GATK.density
    python AncesterAnalysis.py --input Citrus_6strain.vcf.tab --strain $strain
done

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

