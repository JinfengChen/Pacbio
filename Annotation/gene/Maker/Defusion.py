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
python Defusion.py --input ClementinaChr1.all

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


def get_gene_line(prefix, predictor):
    gff = '%s.%s.gff' %(prefix, predictor)
    gff_gene = '%s.%s.gene.gff' %(prefix, predictor)
    ofile = open(gff_gene, 'w')
    with open (gff, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if (unit[2]=="match"):
                    print >> ofile, line
    ofile.close()


'''
summary gene overlap
'''
def summary_overlap(overlap, pd1, pd2):
    test = defaultdict(lambda : int())
    pd_score = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    with open (overlap, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                gene = '%s:%s:%s:%s' %(pd1, unit[0], unit[3], unit[4])
                test[gene] += 1
    for gene in test.keys():
        if test[gene] > 1:
            pd_score[pd1][gene][pd2] = test[gene]
    return pd_score

'''
overlap gff with other predicted gff and clean
'''
def overlap_gff(prefix, pd_list):
    pd_score = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
    for pd1 in pd_list:
        for pd2 in pd_list:
            if not pd1 == pd2:
                overlap='bedtools intersect -a %s.%s.gene.gff -b %s.%s.gene.gff -wo > %s.%s_%s.overlap' %(prefix, pd1, prefix, pd2, prefix, pd1, pd2)
                os.system(overlap)
                temp_score = summary_overlap('%s.%s_%s.overlap' %(prefix, pd1, pd2), pd1, pd2)
                pd_score.update(temp_score)
   
    for pd1 in pd_list:
        gff          = '%s.%s.gff' %(prefix, pd1)
        gff_gene     = '%s.%s.gene.gff' %(prefix, pd1)
        gff_clean    = '%s.%s.gene_clean.gff' %(prefix, pd1)
        gff_defusion = '%s.%s.defusion.gff' %(prefix, pd1) 
        clean_gene(pd1, gff_gene, gff_clean, pd_score) 
        defuse       = 'bedtools intersect -a %s -b %s > %s' %(gff, gff_clean, gff_defusion)
        os.system(defuse)

'''
clean gff by remove gene overlap two or more gene in two other predictor
'''
def clean_gene(pd, gff, gff_clean, pd_score):
    ofile = open(gff_clean, 'w')
    with open (gff, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                gene = '%s:%s:%s:%s' %(pd, unit[0], unit[3], unit[4])
                if not pd_score[pd].has_key(gene):
                    print >> ofile, line
                else:
                    print len(pd_score[pd][gene].keys())
                    print '%s\t%s' %(pd, pd_score[pd][gene])
                #elif len(pd_score[pd][gene].keys()) < 2:
                #    print >> ofile, line
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


    predictor = ["fgenesh", "augustus", "genemark"]
    for pd in predictor:
        get_gene_line(args.input, pd)

    overlap_gff(args.input, predictor)

if __name__ == '__main__':
    main()

