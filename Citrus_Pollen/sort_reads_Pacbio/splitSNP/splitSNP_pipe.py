#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import pysam
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python splitSNP/splitSNP_pipe.py --input scaffold_1.snp.list --bam citrus_canu1_3_ass.bam

    '''
    print message

#SBATCH --nodes=16
#SBATCH --mem=60G
#SBATCH --time=20:00:00
#SBATCH -p highmem
#SBATCH --workdir=./
#SBATCH --output=splitSNP.stdout
def runjob(script, lines):
    #cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ntasks=12,time=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 30 --lines %s --interval 120 --node 1 --mem 20G --time 100:00:00 --convert no %s' %(lines, script)
    os.system(cmd)

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data

def splitSNPsite(snpfile, prefix):
    os.system('grep "^SNP" -v %s |  split -l 10000 - %s.snp.subfile' %(snpfile, prefix)) 
    snpfiles=glob.glob('%s.snp.subfile*' %(prefix))
    return snpfiles

def calculate_reads(reads_files):
    data = defaultdict(lambda : int())
    for infile in reads_files:
        with open (infile, 'r') as filehd:
            for line in filehd:
                line = line.rstrip()
                if len(line) > 2: 
                    unit = re.split(r'\t',line)
                    data[unit[0]] += 1 
    return data 

#python `pwd`/splitSNP/splitSNP_pacbio.py citrus_canu1_3_ass.bam t1 scaffold_1:47220:C:T
def splitBAM_by_SNP(bam, snpfiles, prefix):
    script = os.path.dirname(os.path.abspath(sys.argv[0]))
    shell = open('%s.sh' %(prefix), 'w')
    for snp in sorted(snpfiles):
        cmd = 'python %s/splitSNP_runner.py --input %s --bam %s --output %s' %(script, snp, bam, snp)
        #os.system(cmd)
        print >> shell, cmd
    shell.close()
    runjob('%s.sh' %(prefix), 1)
    refreads_files = glob.glob('%s/*.snp.subfile*.ref.reads' %(prefix))
    altreads_files = glob.glob('%s/*.snp.subfile*.alt.reads' %(prefix))
    refreads = calculate_reads(refreads_files)
    altreads = calculate_reads(altreads_files)
    reads_all = refreads.keys()
    reads_all.extend(altreads.keys())
    
    ##split reads 
    read_ref_file = open('%s.ref.reads.list' %(prefix), 'w')
    read_alt_file = open('%s.alt.reads.list' %(prefix), 'w')
    read_amb_file = open('%s.amb.reads.list' %(prefix), 'w')
    read_all_sum  = open('%s.all.reads.sum' %(prefix), 'w')
    read_ref_list = []
    read_alt_list = []
    read_amb_list = []
    print >> read_all_sum, "Reads\tread_ref_count\tread_alt_count"
    for read in reads_all:
        read_ref_count = refreads[read] if refreads.has_key(read) else 0
        read_alt_count = altreads[read] if altreads.has_key(read) else 0
        if read_ref_count >= read_alt_count*2:
            print >> read_ref_file, read
            read_ref_list.append(read)
        elif read_alt_count >= read_ref_count*2:
            print >> read_alt_file, read
            read_alt_list.append(read)
        else:
            print >> read_amb_file, read
            read_amb_list.append(read)
        print >> read_all_sum, '%s\t%s\t%s' %(read, read_ref_count, read_alt_count) 
    read_ref_file.close()
    read_alt_file.close()
    read_amb_file.close()
    read_all_sum.close()

    ##generate bam file for splited reads
    summary = [0,0,0,0,0]
    samfile = pysam.AlignmentFile(bam, "rb")
    if not samfile._hasIndex():
        print("BAM file '%s' does not have an index, creating one..." % args.input_bam)
        samfile.close()
        pysam.index(args.input_bam)
        samfile = pysam.AlignmentFile(args.input_bam, "rb")
    bam_files = []
    ref_filename = '%s.ref.reads.bam' %(prefix)
    alt_filename = '%s.alt.reads.bam' %(prefix)
    amb_filename = '%s.amb.reads.bam' %(prefix)
    hom_filename = '%s.hom.reads.bam' %(prefix)
    bam_files = [ref_filename, alt_filename, amb_filename, hom_filename]
    ref_bam = pysam.AlignmentFile(ref_filename, "wb", template=samfile)
    alt_bam = pysam.AlignmentFile(alt_filename, "wb", template=samfile)
    amb_bam = pysam.AlignmentFile(amb_filename, "wb", template=samfile)
    hom_bam = pysam.AlignmentFile(hom_filename, "wb", template=samfile)
    for read in samfile.fetch():
        if read.query_name in read_ref_list:
            ref_bam.write(read)
            summary[0] += 1
        elif read.query_name in read_alt_list:
            alt_bam.write(read)
            summary[1] += 1
        elif read.query_name in read_amb_list:
            amb_bam.write(read)
            summary[2] += 1
        else:
            hom_bam.write(read)
            summary[3] += 1
        summary[4] += 1
    samfile.close()
    ref_bam.close()
    alt_bam.close()
    amb_bam.close()
    hom_bam.close()
    pysam.index(ref_filename)
    pysam.index(alt_filename)
    pysam.index(amb_filename)
    pysam.index(hom_filename)
    split_sum_file = '%s.summary' %(prefix)
    split_sum = open(split_sum_file, 'w')
    print >> split_sum, 'Reference reads: %s' %(summary[0])
    print >> split_sum, 'Alternative reads: %s' %(summary[1])
    print >> split_sum, 'Ambiguous reads: %s' %(summary[2])
    print >> split_sum, 'Homozygous reads: %s' %(summary[3])
    print >> split_sum, 'Total reads: %s' %(summary[4])
    split_sum.close()
  
def convert_sequence(bam_files, prefix): 
    ##bam to fastq and fasta
    bedtools = '/opt/linux/centos/7.x/x86_64/pkgs/bedtools/2.25.0/bin/bedtools'
    seqtk    = '/opt/linux/centos/7.x/x86_64/pkgs/seqtk/281537b/seqtk'
    shell = open('%s.convert.sh' %(prefix), 'w') 
    for bam in bam_files:
        bam = os.path.abspath(bam)
        fastq = re.sub(r'.bam', r'.fastq', bam)
        fasta = re.sub(r'.bam', r'.fasta', bam)
        cmd1 = '%s bamtofastq -i %s -fq %s' %(bedtools, bam, fastq)
        cmd2 = '%s seq -A %s > %s' %(seqtk, fastq, fasta)
        print >> shell, cmd1
        print >> shell, cmd2
    shell.close()
    runjob('%s.convert.sh' %(prefix), 2)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-b', '--bam')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'Pacbio_haplotype_reads'  
 
    prefix = args.output
    outdir = os.path.abspath(args.output)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    #cut SNP file 
    bam = os.path.abspath(args.bam) 
    snpfiles = splitSNPsite(args.input, '%s/%s' %(outdir, prefix))
    #splitBAM
    bam_files = []
    if os.path.exists('%s.summary' %(args.output)):
        bam_files=glob.glob('%s.*.reads.bam' %(args.output)) 
    else:
        bam_files= splitBAM_by_SNP(bam, snpfiles, args.output)
    #convert seq format
    convert_sequence(bam_files, args.output) 

if __name__ == '__main__':
    main()

