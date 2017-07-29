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
python ../../Split_Maker_GFF.py --input ClementinaChr1.all.gff

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    prefix = os.path.splitext(infile)[0]
    ofile_mk = open('%s.maker_gene.gff' %(prefix), 'w')
    ofile_gm = open('%s.genemark.gff' %(prefix), 'w')
    ofile_ag = open('%s.augustus.gff' %(prefix), 'w')
    ofile_sn = open('%s.snap.gff' %(prefix), 'w')
    ofile_fg = open('%s.fgenesh.gff' %(prefix), 'w')
    ofile_pr = open('%s.protein2genome.gff' %(prefix), 'w')
    ofile_es = open('%s.est2genome.gff' %(prefix), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'chr') or line.startswith(r'scaf'):
                unit = re.split(r'\t',line)
                if unit[1] == 'genemark':
                    print >> ofile_gm, line
                elif unit[1] == 'augustus_masked':
                    print >> ofile_ag, line
                elif unit[1] == 'snap_masked':
                    print >> ofile_sn, line
                elif unit[1] == 'fgenesh_masked':
                    print >> ofile_fg, line
                elif unit[1] == 'protein2genome':
                    print >> ofile_pr, line
                elif unit[1] == 'est2genome':
                    print >> ofile_es, line
                elif unit[1] == 'maker':
                    print >> ofile_mk, line
    ofile_mk.close()
    ofile_gm.close()
    ofile_ag.close()
    ofile_sn.close()
    ofile_fg.close()
    ofile_pr.close()
    ofile_es.close()
    return data


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

    readtable(args.input)

if __name__ == '__main__':
    main()

