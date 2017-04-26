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

echo "20kb long reads"
python splitSNP/splitSNP_pipe.py --input scaffold_all.snp.list --bam citrus_20kb_bwa_mem.bam --size 5000 > log 2>&1 &

