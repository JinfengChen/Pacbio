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
python Convert_position_to_one_chromosome_bed.py --chr ../Citrus_diversity/bin/Cclementina_v1.0_scaffolds.chrlen --bed test.bed.hist.cov.bed

Convert bed position from 12 chromosome to 1 chromosome.
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
1       0       100000  0.67    67033
1       50000   150000  1.00    100379
1       100000  200000  1.00    100000
1       150000  250000  1.08    107534
1       200000  300000  1.25    125074
1       250000  350000  1.20    120427
1       300000  400000  1.03    102887
'''
def convert_gff(infile, pos):
    data = defaultdict(lambda : int())
    outfile = re.sub(r'.bed$', r'.1chr.bed', infile)
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'scaffold'):
                unit = re.split(r'\t', line)
                unit[1] = str(pos[unit[0]] + int(unit[1]))
                unit[2] = str(pos[unit[0]] + int(unit[2]))
                unit[0] = 'scaffold_1'
                print >> ofile, '\t'.join(unit)
            else:
                print >> ofile, line
    ofile.close()
    return outfile

'''
scaffold_1	28940638
scaffold_2	36376123
scaffold_3	51050279
scaffold_4	25649528
scaffold_5	43300495
scaffold_6	25613772
'''
def read_chr(infile, bed):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'scaffold'): 
                unit = re.split(r'\t',line)
                chrs = int(re.sub(r'scaffold_', r'',unit[0]))
                data[chrs] = int(unit[1])
                #print chrs, data[chrs]
    last = 0
    pos  = defaultdict(lambda : int())
    pos0  = defaultdict(lambda : int())
    pos0[1] = 0
    newchrlen = re.sub(r'.bed$', r'.1chr.chrlen',bed)
    ofile = open(newchrlen, 'w')
    for c in sorted(data.keys(), key=int):
        last  += data[c]
        pos[c] = last
        pos0['scaffold_%s' %(c+1)]= last
        print >> ofile, 'scaffold_%s\t%s' %(1, pos[c])
    ofile.close()
    return pos0

def plot_cov(bed):
    prefix = re.sub(r'.bed$', r'.1chr', bed)
    R_cmd='''
pdf("%s.pdf", width = 10, height = 3)
par(mfrow=c(1, 1))
par(mar=c(5,4,4,2))
x <- read.table("%s.bed")
y <- read.table("%s.chrlen")
plot(x[,2]/1000000, x[,4], col='black', type='l', cex=0.2, main="", xlab="Chromosomal Position (Mb)",  ylab="Read Depth", xlim=c(0, 290), ylim=c(0, 200), frame.plot=FALSE)
for (i in y[,2]){
#abline(v=i/1000000, col='black', lty=3)
segments(i/1000000, 0, i/1000000, 220, col='gray', lty=2, xpd=TRUE)
}
''' %(prefix, prefix, prefix)
    ofile = open('%s.R' %(prefix), 'w')
    print >> ofile, R_cmd
    ofile.close()
    os.system('cat %s.R | R --slave' %(prefix))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chr')
    parser.add_argument('-b', '--bed')
    args = parser.parse_args()
    try:
        len(args.chr) > 0
    except:
        usage()
        sys.exit(2)

    pos = read_chr(args.chr, args.bed)
    outfile = convert_gff(args.bed, pos)
    plot_cov(args.bed)
 
if __name__ == '__main__':
    main()

