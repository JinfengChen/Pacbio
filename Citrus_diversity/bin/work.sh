echo "chromosome length"
head -n 9 ~/BigData/00.RD/Assembly/Pacbio/Reference/Cclementina_JGI_v1.0/Cclementina_v1.0_scaffolds.chrlen > Cclementina_v1.0_scaffolds.chrlen
python Convert_position_to_one_chromosome_tab.py --chr Cclementina_v1.0_scaffolds.chrlen --tab Cclementina_v1.0_scaffolds.chrlen > Cclementina_v1.0_scaffolds.chrlen.1chr

echo "vcf to tab to heterozygous rate"
module load vcftools
zcat ../input/SNP/SWO.vcf.gz| vcf-to-tab > SWO.vcf.tab
python HeterozygousDensity.py --input SWO.vcf.tab > SWO.vcf.tab.het.density &
python Convert_position_to_one_chromosome_tab.py --chr Cclementina_v1.0_scaffolds.chrlen --tab SWO.vcf.tab.het.density > SWO.vcf.tab.het.1chr.density
python plot_chr_density.py --input SWO.vcf.tab.het.1chr.density

#run for all strain by run.sh
qsub run.sh

echo "Ancestral analysis"
cat ../input/SNP/Citrus_4strain.vcf | vcf-to-tab > Citrus_4strain.vcf.tab
python HeterozygousDensity.py --input Citrus_4strain.vcf.tab --strain Citrus > FAIRCHILD.vcf.tab.het.GATK.density &
python Convert_position_to_one_chromosome_tab.py --chr Cclementina_v1.0_scaffolds.chrlen --tab FAIRCHILD.vcf.tab.het.GATK.density > FAIRCHILD.vcf.tab.het.1chr.GATK.density
python plot_chr_density.py --input FAIRCHILD.vcf.tab.het.1chr.GATK.density
python AncesterAnalysis.py --input Citrus_4strain.vcf.tab --strain Citrus &

