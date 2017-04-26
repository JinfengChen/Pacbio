splitSNP
========

*splitSNP* is a small utility to extract reads from a BAM file that cover a specified SNP and output the reference and alternate allele containing reads to separate BAM files.

I created *splitSNP* for the analysis of allele-specific methylation from bisulfite amplicon experiments sequenced on the MiSeq.

Usage
---------
    splitSNP.py [-h] [--max_depth MAX_DEPTH] input_bam output_prefix SNP


* input_bam - Input BAM file
* output_prefix - Prefix for output files
* SNP - SNP position, reference and alternate allele of interest in chr:position:ref:alt format, eg chr21:11106932:A:G. For deletion analysis the ref should be 'D' and alt be the size of the deletion in basepairs, eg chr11:67351213:D:64 .Coordinates are 1-based.
* --pair_distance PAIR_DISTANCE The distance in basepairs to search up and downstream from the specified SNP/deletion for the pair of overlapping reads (default is 500)
* --max_depth MAX_DEPTH -  Maximum number of reads to process at the specified SNP position (default is 1000000)

SNP usage
---------
This example uses *splitSNP* to extract and separate the **A** (reference allele) and **G** (alternate allele) containing reads spanning the SNP *rs11184058* which located at chr21, position 11,106,932 in hg19 (obtained from UCSC genome browser).

    splitSNP.py MiSeq_20.bam MiSeq_20.rs11184058 chr21:11106932:A:G

<!-- separate input and output -->

    BAM file 'MiSeq_20.bam' does not have an index, creating one...
    Processing 118528 reads covering SNP position chr21:11106932 in MiSeq_20.bam
    20 reads discarded for being ambiguous
    187646 read segments with the reference allele 'A' written to MiSeq_20.rs11184058.ref.A.bam
    11026 read segments with the alternate allele 'G' written to MiSeq_20.rs11184058.alt.G.bam
    Discarded 237 'C' alleles
    Discarded 42 'T' alleles
    Discarded 35 'N' alleles

NOTE: 118,528 reads overlap the SNP *rs11184058*, however 198,672 reads are written out to the two output files as read pairs within the basepair distance specifed by the *--pair_distance* parameter (default 500) are also outputted.

Deletion usage
--------------
This example uses *splitSNP* to extract reads from the human *GSTP1* locus, separating reads matching the reference sequence from those containing an engineered 64 basepair deletion spanning chr11:67,351,213-67,351,276 (hg19).



    splitSNP.py MiSeq_1_S1.bam MiSeq_1_S1.GSTP1 chr11:67351213:D:64

<!-- separate input and output -->

    13213 read segments with the reference sequence written to GSTP1_Split/MiSeq_1_S1/MiSeq_1_S1.ref.bam
    4877 read segments with the 64bp deletion written to GSTP1_Split/MiSeq_1_S1/MiSeq_1_S1.del.64bp.bam


TODOs
-----
* There's nothing here!

CHANGELOG
---------
* 2014-01-12 - Added *--pair_distance* parameter so that read pairs that do not overlap the SNP/deletion are also output
* 2014-01-8 - Added deletion mode
* 2014-01-8 - Created initial version
