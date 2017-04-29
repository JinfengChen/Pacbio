echo "testing split SNP pacbio pipeline"
python splitSNP/splitSNP_pipe.py --input scaffold_1.snp.list --bam citrus_canu1_3_ass.bam > log 2>&1 &
echo "preads split scaffold_1"
python splitSNP/splitSNP_pipe.py --input scaffold_1.snp.list --bam citrus_preads_bwa_mem_scaffold1.bam --size 500 > log 2>&1 &
echo "preads split scaffold_all"
python splitSNP/splitSNP_pipe.py --input scaffold_all.snp.list --bam citrus_preads_bwa_mem.bam --size 5000 > log 2>&1 &
echo "rawreads split scaffold_all"
python splitSNP/splitSNP_pipe.py --input scaffold_all.snp.list --bam citrus_rawreads_bwa_mem.bam --size 5000 > log 2>&1 &

echo "extract reads from single fatsa file"
python splitSNP/splitSNP_extract_reads_Pacbio.py --input Pacbio_haplotype_reads --fasta citrus.pacbio.fasta
echo "extract reads from list of fatsa files"
python splitSNP/splitSNP_extract_reads_Pacbio.py --input Pacbio_haplotype_reads --fasta_list input.fofn

echo "20kb long reads, need to use pysam < 0.8.3"
#The _hasIndex() call in the __init__.py script is broken due to changes in the pysam.csamfile.Samfile object in pysam version >0.8.3 which has renamed this function to has_index().
ln -s /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Testing_Data/Citrus/SMARTcell_20kb/citrus.pacbio_20kb.fasta ./
perl ~/BigData/software/bin/fastaDeal.pl --attr id citrus.pacbio_20kb.fasta > citrus.pacbio_20kb.fasta.id
grep "20kb" /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON_v3.0/virtualenv/FALCON_v3.0/Citrus/Pacbio_raw_plus_20kb/input.fofn > input.20kb.fofn
python splitSNP/splitSNP_pipe.py --input scaffold_all.snp.list --bam citrus_20kb_bwa_mem.bam --size 5000 --output Pacbio_20kb_haplotype_reads
python splitSNP/splitSNP_extract_reads_Pacbio.py --input Pacbio_20kb_haplotype_reads --fasta citrus.pacbio_20kb.fasta --fastaid citrus.pacbio_20kb.fasta.id
python splitSNP/splitSNP_extract_reads_Pacbio.py --input Pacbio_20kb_haplotype_reads --fasta_list input.20kb.fofn --fastaid citrus.pacbio_20kb.fasta.id
sbatch splitSNP.sh
