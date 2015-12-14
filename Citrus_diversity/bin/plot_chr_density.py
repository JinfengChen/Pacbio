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

def plot_density(prefix, density, chromosome):
    R='''
pdf("%s.pdf", width = 10, height = 4)
x <- read.table("%s", header=TRUE)
y <- read.table("%s")
plot(x[,2]/1000000, x[,5], col='cadetblue', type='l', main="", xlab="Chromosomal Position (Mb)",  ylab="Heterozygous site/kb", xlim=c(0, 290), ylim=c(0, 30), frame.plot=FALSE)
for (i in y[,2]){
abline(v=i/1000000, col='black', lty=3)
}
dev.off()
''' %(prefix, density, chromosome)
    ofile = open('%s.R' %(prefix), 'w')
    print >> ofile, R
    ofile.close()
    os.system('cat %s.R | R --slave' %(prefix))

def plot_distri(prefix, density):
    R='''
pdf("%s.distri.pdf", width = 10, height = 4)
x <- read.table("%s", header=TRUE)
plot(density(x[,5]), xlab="Heterozygous site/kb", ylab="Frequency")
dev.off()
''' %(prefix, density)
    ofile = open('%s.distri.R' %(prefix), 'w')
    print >> ofile, R
    ofile.close()
    os.system('cat %s.distri.R | R --slave' %(prefix))
    

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
    prefix = os.path.splitext((os.path.basename(args.input)))[0]
    #print prefix
    plot_density(prefix, args.input, 'Cclementina_v1.0_scaffolds.chrlen.1chr')
    plot_distri(prefix, args.input)

if __name__ == '__main__':
    main()

