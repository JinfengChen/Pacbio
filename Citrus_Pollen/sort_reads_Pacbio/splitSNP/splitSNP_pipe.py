#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
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
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 30 --lines %s --interval 120 --node 2 --mem 20G --time 100:00:00 --convert no %s' %(lines, script)
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
    
    bam = os.path.abspath(args.bam) 
    snpfiles = splitSNPsite(args.input, '%s/%s' %(outdir, prefix))
    splitBAM_by_SNP(bam, snpfiles, args.output)
     

if __name__ == '__main__':
    main()

