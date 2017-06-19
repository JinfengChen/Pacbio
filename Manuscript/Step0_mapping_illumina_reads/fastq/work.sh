python SRA_down.py --input sra.list > log 2>&1 &
python SRA_merge.py --input ./Citrus_RNAseq > log 2>&1 &
mv Citrus_RNAseq/SRR1023*.gz ./
bash rename.sh
mv ~/BigData/00.RD/Assembly/Pacbio/Citrus_diversity/input/fastq/30X/p00.CLM_* ./
sbatch pgz.sh

