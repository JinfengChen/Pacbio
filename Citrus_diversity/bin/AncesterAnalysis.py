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
#CHROM  POS     REF     CLM     LAP     WMM
scaffold_4      830     A       A/A     ./.     A/C
scaffold_4      1234    T       T/T     ./.     T/G
scaffold_4      1531    T       T/T     ./.     T/C
scaffold_4      1861    G       G/G     ./.     G/T
scaffold_4      1862    TGGAAGG TGGAAGG/TGGAAGG ./.     TGGAAGG/T
'''
def read_tab(infile, strain):
    data = defaultdict(lambda : defaultdict(lambda : list()))
    columen = 3
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'scaffold'): 
                unit = re.split(r'\t',line)
                if len(unit[2]) == 1 and len(unit[columen]) == 3:
                    #print '%s\t%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[3], unit[3][0], unit[3][-1])
                    data[unit[0]][unit[1]] = [unit[columen][0], unit[columen][-1]]
            elif line.startswith('#CHROM'):
                unit = re.split(r'\t',line)
                for i in range(len(unit)):
                    if unit[i] == strain:
                        columen = i
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

def window_analysis(aSNP, SNP, prefix):
    win  = 2000
    step = 2000
    ofile = open('%s.ancester.bed' %(prefix), 'w')
    #create scaffold_1 to scaffold 9 to loop
    for i in range(1, 10):
    #for i in [4]:
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
                    #print '%s\t%s\tC_maxima:%s\tC_reticulata:%s\t%s\t%s' %(scaf, aSNP_scaf_win[rank][1], aSNP_scaf_win[rank][2], aSNP_scaf_win[rank][3], SNP[scaf][aSNP_scaf_win[rank][1]][0], SNP[scaf][aSNP_scaf_win[rank][1]][1])
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
            #print '%s\t%s\t%s' %(ancester[0], ancester[1], ancester[2])
            ancester_state = 'unknown'
            ancester_snp   = 0
            #C_maxima
            if ancester[0] >= ancester[1] and ancester[0] >= len(aSNP_scaf_win)/2:
                ancester_state = 'C_maxima'
                ancester_snp   = '%s/%s/%s' %(ancester[0], ancester[1], ancester[2])
                #ancester_snp   = ancester[0]
            #C_reticulata
            elif ancester[1] >= ancester[0] and ancester[1] >= len(aSNP_scaf_win)/2:
                ancester_state = 'C_reticulata'
                ancester_snp   = '%s/%s/%s' %(ancester[0], ancester[1], ancester[2])
                #ancester_snp   = ancester[1]
            #C_maxima/C_reticulata
            elif ancester[2] >= ancester[0] and ancester[2] >= ancester[1] and ancester[2] >= len(aSNP_scaf_win)/2:
                ancester_state = 'C_maxima/C_reticulata'
                ancester_snp   = '%s/%s/%s' %(ancester[0], ancester[1], ancester[2])
                #ancester_snp   = ancester[2]
            else:
                ancester_state = 'unknown'
                ancester_snp   = '%s/%s/%s' %(ancester[0], ancester[1], ancester[2])
            print >> ofile, '%s\t%s\t%s\t%s\t%s\t%s' %(scaf, aSNP_scaf_win[0][1], aSNP_scaf_win[-1][1], ancester_state, ancester_snp, '+')
    ofile.close()

def merge_overlap_win(prefix):
    raw_bed  = '%s.ancester.bed' %(prefix)
    uniq_bed = '%s.ancester.uniq.bed' %(prefix)
    unique   = []
    with open (raw_bed, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            unit = re.split(r'\t', line)
            if len(unique) == 0:
                unique.append(unit)
            else:
                # current window overlap with last window and with same genotype
                #print '%s\t%s\t%s\t%s' %(unit[0], unique[-1][0], unit[1], unique[-1][2])
                if unit[0] == unique[-1][0] and int(unit[1]) <= int(unique[-1][2]) and unit[3] == unique[-1][3]:
                    unique[-1][2] = int(unit[2])
                else:
                    unique.append(unit)
    ofile = open (uniq_bed, 'w')
    for win in unique:
        if win[3] == 'C_reticulata':
            win[4] = 'orange'
        elif win[3] == 'C_maxima':
            win[4] = 'green'
        elif win[3] == 'C_maxima/C_reticulata':
            win[4] = 'blue'
        else:
            win[4] = 'gray'
        print >> ofile, '\t'.join(map(str, win)) 
    ofile.close() 

'''
scaffold_1      9791    28591350        C_reticulata    1906    +
scaffold_1      28398652        28882450        C_maxima/C_reticulata   1373    +
scaffold_1      28750435        28933880        unknown 0       +
'''
def convert_gff(infile, pos):
    data = defaultdict(lambda : int())
    ofile = open(re.sub(r'.bed', r'.1chr.bed', infile), 'w')
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

#Chr01	43270923	16701176	17133774
#Chr10	23207287	8100966	8178267
def read_chr(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'scaffold'): 
                unit = re.split(r'\t',line)
                chrs = int(re.sub(r'scaffold_', r'',unit[0]))
                data[chrs] = int(unit[1])
    last = 0
    pos  = defaultdict(lambda : int())
    for c in sorted(data.keys(), key=int):
        pos['scaffold_%s' %(c)] = last
        last  += data[c]
        #print 'Chr%s\t%s' %(c, pos[c])
    return pos

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
  
    if not args.output:
        args.output = args.strain
        if args.output == 'Citrus':
            args.output = 'FAIRCHILD'

    if not os.path.exists('%s.ancester.bed' %(args.output)):
        aSNP = ancester_SNP('../input/nbt.2906-S2.txt')
        SNP  = read_tab(args.input, args.strain) 
        window_analysis(aSNP, SNP, args.output)
    merge_overlap_win(args.output)
    pos = read_chr('Cclementina_v1.0_scaffolds.chrlen')
    convert_gff('%s.ancester.uniq.bed' %(args.output), pos)
 
if __name__ == '__main__':
    main()

