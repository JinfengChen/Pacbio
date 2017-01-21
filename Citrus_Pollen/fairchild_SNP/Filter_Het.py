#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
import scipy.stats
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python Filter_Het.py --input dbSNP.clean.vcf

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

###source=SelectVariants
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Citrus
#scaffold_1	33023	.	T	C	93645.4	PASS	.	GT:AD:DP:GQ:PL	1/1:0,30:30:90:1043,90,0
def Filter_Het(infile):
    data = defaultdict(str)
    HetSNP_vcf = re.sub(r'.vcf', r'.Het.vcf', infile)
    ofile = open(HetSNP_vcf, 'w')
    print HetSNP_vcf
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                info = re.split(r':', unit[9])
                if info[0] == '0/1' and info[2] >= 4: #het and have more than 4 reads covered
                    #print info[0], line
                    alleles = re.split(r',', info[1])
                    alleles.sort(key=int)
                    #print alleles[0], info[2]
                    pvalue  = scipy.stats.binom_test(alleles[0],n=info[2],p=0.5,alternative='less')
                    if pvalue > 0.05:
                        print >> ofile, line, 'good'
                    else:
                        print >> ofile, line, 'bad'
            else:
                #print >> ofile, line
                pass
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    Filter_Het(args.input)

if __name__ == '__main__':
    main()

