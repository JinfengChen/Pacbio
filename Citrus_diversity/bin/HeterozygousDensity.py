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
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
#CHROM  POS     REF     CLM     LAP     WMM
scaffold_4      830     A       A/A     ./.     A/C
scaffold_4      1234    T       T/T     ./.     T/G
scaffold_4      1531    T       T/T     ./.     T/C
scaffold_4      1861    G       G/G     ./.     G/T
scaffold_4      1862    TGGAAGG TGGAAGG/TGGAAGG ./.     TGGAAGG/T
'''
def read_tab(infile, strain):
    data = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    window=200000
    columen = 3
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'scaffold'): 
                unit = re.split(r'\t',line)
                #only count SNP not indel
                if len(unit[2]) == 1 and len(unit[columen]) == 3:
                    index = int(unit[1])//int(window)
                    midpoint = int(index)*window + window/2 
                    data[unit[0]][midpoint][0] += 1
                    #only if heterozygous SNP
                    if not unit[columen][0] == unit[columen][-1]:
                        data[unit[0]][midpoint][1] += 1
            elif line.startswith('#CHROM'):
                unit = re.split(r'\t',line)
                for i in range(len(unit)):
                    #strain not specific, use default value for columen
                    if strain == 'NA':
                        continue
                    #strain specific, use strain value
                    elif unit[i] == strain:
                        columen = i
    for num in range(1, 10):
        scaf = 'scaffold_%s' %(num)
        #print scaf
        for pos in sorted(data[scaf].keys(), key=int):
            het_rate = float(data[scaf][pos][1])/window*1000
            print '%s\t%s\t%s\t%s\t%s' %(scaf, pos, data[scaf][pos][0], data[scaf][pos][1], het_rate)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-s', '--strain')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.strain:
        strain = 'NA'

    read_tab(args.input, args.strain)

if __name__ == '__main__':
    main()

