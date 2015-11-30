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
python AncesterAnalysis.py --input ../input/nbt.2906-S2.txt

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

'''
#CHROM  POS     REF     CLM
scaffold_1      54      CC      CA/CA
scaffold_1      75      T       T/C
scaffold_1      164     A       A/G
scaffold_1      302     C       C/G
scaffold_1      330     T       T/G
'''
def read_tab(infile):
    data = defaultdict(lambda : defaultdict(lambda : list()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'scaffold'): 
                unit = re.split(r'\t',line)
                if len(unit[2]) == 1 and len(unit[3]) == 3:
                    #print '%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[3], unit[3][0], unit[3][-1])
                    data[unit[0]][unit[1]] = [unit[3][0], unit[3][-1]]
    return data



##chromosome  position  C_maxima   C_reticulata
#scaffold_1      9791    G       A
#scaffold_1      38000   A       C
#scaffold_1      38040   A       G
#scaffold_1      38042   G       A
#scaffold_1      38057   C       T
def ancester_SNP(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]].append(unit)
    return data

def window_analysis(aSNP, SNP):
    win  = 2000
    step = 1000
    for i in range(1, 10):
        scaf = 'scaffold_%s' %(i)
        #print scaf
        aSNP_scaf = aSNP[scaf]
        cycles = len(aSNP_scaf)/step
        #sliding window
        for w in range(cycles):
            start_w = w*step
            end_w   = w*step + win
            aSNP_scaf_win = aSNP_scaf[start_w:end_w]
            #print w, len(aSNP_scaf_win)
            #print aSNP_scaf_win[0], aSNP_scaf_win[0][1]
            #print aSNP_scaf_win[-1], aSNP_scaf_win[-1][1]
            #C_maxima, C_reticulata, C_maxima/C_reticulata
            ancester = [0, 0, 0]
            #analyze ancestral SNP in each window
            for rank in range(len(aSNP_scaf_win)):
                #analyze only these SNP overlapped with diagnostic SNP
                if SNP[scaf].has_key(aSNP_scaf_win[rank][1]):
                    print '%s\t%s\tC_maxima:%s\tC_reticulata:%s\t%s\t%s' %(scaf, aSNP_scaf_win[rank][1], aSNP_scaf_win[rank][2], aSNP_scaf_win[rank][3], SNP[scaf][aSNP_scaf_win[rank][1]][0], SNP[scaf][aSNP_scaf_win[rank][1]][1])
                    #homozygous
                    if SNP[scaf][aSNP_scaf_win[rank][1]][0] == SNP[scaf][aSNP_scaf_win[rank][1]][1]:
                        #C_maxima
                        if SNP[scaf][aSNP_scaf_win[rank][1]][0] == aSNP_scaf_win[rank][2]:
                            ancester[0] += 1
                        #C_reticulata
                        elif SNP[scaf][aSNP_scaf_win[rank][1]][0] == aSNP_scaf_win[rank][3]:
                            ancester[1] += 1
                    #heterzygous
                    else:
                        #C_maxima/C_reticulata
                        if SNP[scaf][aSNP_scaf_win[rank][1]][0] == aSNP_scaf_win[rank][2] and SNP[scaf][aSNP_scaf_win[rank][1]][1] == aSNP_scaf_win[rank][3]:
                            ancester[2] += 1
                        #C_maxima/C_reticulata
                        elif SNP[scaf][aSNP_scaf_win[rank][1]][0] == aSNP_scaf_win[rank][3] and SNP[scaf][aSNP_scaf_win[rank][1]][1] == aSNP_scaf_win[rank][2]:
                            ancester[2] += 1
            print '%s\t%s\t%s' %(ancester[0], ancester[1], ancester[2])
            ancester_state = 'unknown'
            ancester_snp   = 0
            #C_maxima
            if ancester[0] >= ancester[1] and ancester[0] >= len(aSNP_scaf_win)/2:
                ancester_state = 'C_maxima'
                ancester_snp   = ancester[0]
            #C_reticulata
            elif ancester[1] >= ancester[0] and ancester[1] >= len(aSNP_scaf_win)/2:
                ancester_state = 'C_reticulata'
                ancester_snp   = ancester[1]
            #C_maxima/C_reticulata
            elif ancester[2] >= ancester[0] and ancester[2] >= ancester[1] and ancester[2] >= len(aSNP_scaf_win)/2:
                ancester_state = 'C_maxima/C_reticulata'
                ancester_snp   = ancester[2]
            print '%s\t%s\t%s\t%s\t%s\t%s' %(scaf, aSNP_scaf_win[0][1], aSNP_scaf_win[-1][1], ancester_state, ancester_snp, '+')

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

    aSNP = ancester_SNP('../input/nbt.2906-S2.txt')
    SNP  = read_tab(args.input) 
    #aSNP = '' 
    window_analysis(aSNP, SNP)

if __name__ == '__main__':
    main()

