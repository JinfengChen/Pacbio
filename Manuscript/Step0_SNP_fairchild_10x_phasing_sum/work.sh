wget https://gist.githubusercontent.com/slowkow/6215557/raw/252541a75b531d974bf5a9eae1b0d0a85cce225d/VCF.py
cp ../Step0_SNP_fairchild_10x_phasing/fairchild_wgs_phasing/outs/phased_variants.vcf.gz ./

python Phased_SNP_sum.py --vcf phased_variants.vcf.gz --haplotype pollen_haplotype/Haplotype.FCM.txt > log 2>&1 &
awk '$2>10' phased_variants.vcf.gz.block_sum.txt | awk '{print $6/$2, $7/$2}'| less -S
#compare with pollen
awk '$2>=10' phased_variants.vcf.gz.block_sum.txt | awk '{if ($6+$7 > 0){ if($6>$7){ print $6/($6+$7) }else{ print $7/($6+$7) }  }}' > phased_variants.vcf.gz.block_sum.txt.max_rate.txt

#183 Mb block are larger than 1Mb.
awk '$2>=10' phased_variants.vcf.gz.block_sum.txt | awk '$3>1000000'| cut -f3| perl ~/BigData/software/bin/numberStat.pl
#218 Mb block in total
awk '$2>=10' phased_variants.vcf.gz.block_sum.txt | cut -f3| perl ~/BigData/software/bin/numberStat.pl
#8 Mb block lower than 0.98
awk '$2< 0.98' phased_variants.vcf.gz.block_sum.txt.max_rate.txt | perl ~/BigData/software/bin/numberStat.pl

