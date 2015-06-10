#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
perl CompareReadLength.py --fasta ecoli/preads4falcon.fasta --list ecoli/input.fofn --project Pacbio_FALCON

Generate read length file for fasta file or list of fasta files. Compare read length distribution when both files are given

    '''
    print message


def fasta_len(fastafile, outfile):
    fastaid = defaultdict(str)
    ofile = open(outfile, 'a')
    for record in SeqIO.parse(fastafile,"fasta"):
        #fastaid[record.id] = str(len(str(record.seq)))
        #print >> ofile, '%s\t%s' %(str(record.id), str(len(str(record.seq))))
        print >> ofile, str(len(str(record.seq)))
    ofile.close()
    #return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = 1
    return data

def Draw_Length(title):
    Rcmd='''pdf("%s.Length_Compare.pdf")
corrected <- read.table("%s.corrected.length")
raw    <- read.table("%s.raw.length")
cor_xh <- hist(corrected[,1], plot=FALSE)
raw_xh <- hist(raw[,1], plot=FALSE)
xmax <- max(max(raw_xh$breaks), max(cor_xh$breaks))
ymax <- max(max(raw_xh$counts), max(cor_xh$counts))
atxvalue <- seq (0, xmax, by=10000)
atyvalue <- seq (0, ymax, by=4000)
atx <- c(atxvalue, xmax+20)
aty <- c(atyvalue, ymax+10)
plot(raw_xh$mids, raw_xh$counts, type='b', pch=18, col='blue', xlab="Read Length (bp)", ylab="#Read", xlim=c(0, xmax), ylim=c(0, ymax), axes=FALSE)
lines(cor_xh$mids, cor_xh$counts, type='b', pch=10, col='orange')
axis(1,at=atx,labels=atx)
axis(2,at=aty,labels=aty)
legend("topright", bty="n",lty=c(1, 1), pch=c(18, 10), cex=1.2, c("Raw", "Corrected"), col=c("Blue" ,"Orange"))
dev.off()
''' %(title, title, title)
    ofile = open('%s.Length_Compare.R' %(title), 'w')
    print >> ofile, Rcmd
    ofile.close()
    os.system('cat %s.Length_Compare.R | R --slave' %(title))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta')
    parser.add_argument('-l', '--list')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    if not args.project:
        args.output = 'Pacbio'

    ##corrected read length, normally in one fasta file
    if args.fasta:
        corrected_read_length = '%s.corrected.length' %(args.project)
        if not os.path.isfile(corrected_read_length):
            fasta_len(args.fasta, corrected_read_length)

    ##raw read length, from list of fasta file
    if args.list:
        raw_read_length = '%s.raw.length' %(args.project) 
        raw_reads = readtable(args.list)
        if not os.path.isfile(raw_read_length):
            for readfile in raw_reads:
                fasta_len(readfile, raw_read_length)
 
    if args.fasta and args.list:
        Draw_Length(args.project)

if __name__ == '__main__':
    main()

