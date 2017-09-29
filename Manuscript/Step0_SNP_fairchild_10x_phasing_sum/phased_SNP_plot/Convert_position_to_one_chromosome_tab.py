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
python Convert_position_to_one_chromosome.py --chr Clementine.chr.inf --tab SWO.vcf.tab.het.density 

Convert tab position from 12 chromosome to 1 chromosome.
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#scaffold_1      100000  2874    2186    10.93
#scaffold_1      300000  1439    775     3.875
#scaffold_1      500000  1192    654     3.27
#scaffold_1      700000  1517    1137    5.685
def convert_gff(infile, pos):
    data = defaultdict(lambda : int())
    #ofile = open(re.sub(r'.density', r'.1chr.density', infile), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'chr'):
                unit = re.split(r'\t', line)
                unit[1] = str(pos[unit[0]] + int(unit[1]))
                unit[2] = str(pos[unit[0]] + int(unit[2]))
                unit[0] = 'chr1'
                print '\t'.join(unit)
            else:
                print line
    #ofile.close()

#Chr01	43270923	16701176	17133774
#Chr10	23207287	8100966	8178267
def read_chr(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'chr'): 
                unit = re.split(r'\t',line)
                chrs = int(re.sub(r'chr', r'',unit[0]))
                data[chrs] = int(unit[1])
    last = 0
    pos  = defaultdict(lambda : int())
    for c in sorted(data.keys(), key=int):
        pos['chr%s' %(c)] = last
        last  += data[c]
        #print 'Chr%s\t%s' %(c, pos[c])
    return pos


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chr')
    parser.add_argument('-t', '--tab')
    args = parser.parse_args()
    try:
        len(args.chr) > 0
    except:
        usage()
        sys.exit(2)

    pos = read_chr(args.chr)
    convert_gff(args.tab, pos)

if __name__ == '__main__':
    main()

