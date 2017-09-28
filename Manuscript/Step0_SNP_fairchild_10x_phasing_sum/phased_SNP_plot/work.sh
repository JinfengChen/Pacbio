#python Prepare_pollen.py --input Haplotype.FCM.txt
#head -n 9 ~/BigData/00.RD/Assembly/Pacbio/Reference/Fairchild/Fairchild_v1.fasta.len > Fairchild_v1.fasta.chrlen
#python Convert_position_to_one_chromosome_tab.py --chr Fairchild_v1.fasta.chrlen --tab Fairchild_v1.fasta.chrlen > Fairchild_v1.fasta.1chr.chrlen
python Convert_position_to_one_chromosome_tab.py --chr Fairchild_v1.fasta.chrlen --tab GT.hap1.txt > GT.hap1.1chr.txt
python Convert_position_to_one_chromosome_tab.py --chr Fairchild_v1.fasta.chrlen --tab GT.hap2.txt > GT.hap2.1chr.txt
python Convert_position_to_one_chromosome_tab.py --chr Fairchild_v1.fasta.chrlen --tab GT.Pollen41.txt > GT.Pollen41.1chr.txt
python Convert_position_to_one_chromosome_tab.py --chr Fairchild_v1.fasta.chrlen --tab GT.Pollen42.txt > GT.Pollen42.1chr.txt
python Convert_position_to_one_chromosome_tab.py --chr Fairchild_v1.fasta.chrlen --tab GT.Pollen43.txt > GT.Pollen43.1chr.txt
python Convert_position_to_one_chromosome_tab.py --chr Fairchild_v1.fasta.chrlen --tab GT.Pollen44.txt > GT.Pollen44.1chr.txt
python Convert_position_to_one_chromosome_tab.py --chr Fairchild_v1.fasta.chrlen --tab GT.Pollen45.txt > GT.Pollen45.1chr.txt
cat FAIRCHILD.1chr.GATK.plot_genotype.R | R --slave

