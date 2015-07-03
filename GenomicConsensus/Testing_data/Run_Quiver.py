#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob

def usage():
    test="name"
    message='''
python Run_Quiver.py --cmddir test_cmd_sort_split --ref yeast_ass_line50_clean.fa --project yeast_quiver

#0. modify -j4 and -p in quiver command: -j cpu, -p chemical
#1. ref fasta need to 50 or some other base per line, not single long line per sequence
#2. ref fasta name need to short and do not have space
#3. do not use one SMRTcell to merge and sort
#4. --mapQvThreshold=0 in denovo assembly consensus to deal with low quality read mapping in repeat regions

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 40 --lines %s --interval 120 --resource nodes=1:ppn=4,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

def writefile(outfile, lines):
    ofile = open(outfile, 'w')
    print >> ofile, lines
    ofile.close()


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cmpdir')
    parser.add_argument('-r', '--ref')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.cmpdir) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'Quiver_outdir' 
  
    if not os.path.exists(args.project):
        os.mkdir(args.project)

    #run quiver on each contig
    cmp_files = glob.glob('%s/*.cmp.h5' %(args.cmpdir))
    cmd_lines = []
    for h5 in sorted(cmp_files):
        prefix = os.path.split(h5)[1]
        prefix = re.sub(r'\.cmp\.h5', r'', prefix)
        prefix = '%s/%s' %(args.project, prefix)
        #-p C2.NoQVsModel
        cmd = 'quiver -j4 %s -r %s -o %s.gff -o %s.consensus.fasta -o %s.consensus.fastq' %(os.path.abspath(h5), os.path.abspath(args.ref), os.path.abspath(prefix), os.path.abspath(prefix), os.path.abspath(prefix))
        #print cmd
        cmd_lines.append(cmd)
    writefile('%s.quiver_run.sh' %(args.project), '\n'.join(cmd_lines))
    runjob('%s.quiver_run.sh' %(args.project), 5)

    #merge contig
    merge = []
    merge.append('cat %s/*.consensus.fasta > %s.consensus.fasta' %(os.path.abspath(args.project), os.path.abspath(args.project)))
    merge.append('cat %s/*.consensus.fastq > %s.consensus.fastq' %(os.path.abspath(args.project), os.path.abspath(args.project)))
    merge.append('cat %s/*.gff > %s.gff' %(os.path.abspath(args.project), os.path.abspath(args.project)))
    writefile('%s.quiver_merge.sh' %(args.project), '\n'.join(merge))
    runjob('%s.quiver_merge.sh' %(args.project), 5)

if __name__ == '__main__':
    main()

