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
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines 2 --interval 120 --task 1 --mem 15G --time 100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#ref	alt	pol41	pol42	pol43	pol44	pol45	hap1	hap2
#chr1:557344	C	A	A	A	C	NA	A	C	A
#chr1:557377	T	A	A	A	T	NA	A	T	A
def pollen_haplotype(infile):
    data = defaultdict(lambda : defaultdict(lambda : list()))
    with open (infile, 'r') as filehd:
        pollens = ['Pollen41', 'Pollen42', 'Pollen43', 'Pollen44', 'Pollen45']
        chrs_last = ''
        pos_last  = 0
        hap_chr_last    = defaultdict(lambda : defaultdict(lambda : int()))
        hap_chr_count   = defaultdict(lambda : defaultdict(lambda : int()))
        hap_chr_data    = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : list())))
        hap_chr_hap     = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : int())))
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                chrs, pos = re.split(r':', unit[0])
                if not chrs == chrs_last:
                    #haplotype
                    data['hap1']['start'].append([chrs, pos, 1])
                    data['hap2']['start'].append([chrs, pos, 2])
                    if not chrs_last == '':
                        data['hap1']['end'].append([chrs_last, pos_last])
                        data['hap2']['end'].append([chrs_last, pos_last])
                    #pollen
                    for i in range(3, 8):
                        pol = pollens[i-3]
                        hap = 0
                        if unit[i] == unit[8]:
                            hap = 1
                        elif unit[i] == unit[9]:
                            hap = 2
                        else:
                            hap = 3
                        if hap == 1 or hap == 2:
                            hap_chr_count[pol][chrs] += 1
                            hap_chr_data[pol][chrs][hap_chr_count[pol][chrs]].append(int(pos))
                            hap_chr_hap[pol][chrs][hap_chr_count[pol][chrs]] = hap
                            hap_chr_last[pol][chrs]   = hap
                    chrs_last = chrs
                    pos_last  = pos    
                else:
                    #pollen
                    for i in range(3, 8):
                        pol = pollens[i-3]
                        hap = 0
                        if unit[i] == unit[8]:
                            hap = 1
                        elif unit[i] == unit[9]:
                            hap = 2
                        else:
                            hap = 3

                        if hap == 3 or hap == 0:
                            pass
                        elif hap == hap_chr_last[pol][chrs]:
                            hap_chr_data[pol][chrs][hap_chr_count[pol][chrs]].append(int(pos))
                            hap_chr_hap[pol][chrs][hap_chr_count[pol][chrs]] = hap
                        else:
                            hap_chr_count[pol][chrs] += 1
                            hap_chr_data[pol][chrs][hap_chr_count[pol][chrs]].append(int(pos))
                            hap_chr_hap[pol][chrs][hap_chr_count[pol][chrs]] = hap
                            hap_chr_last[pol][chrs]   = hap 
                        #data[pol]['end'].append([chrs_last, pos_last])
                        #data[pol]['start'].append([chrs, pos, hap])
                    chrs_last = chrs
                    pos_last  = pos
        data['hap1']['end'].append([chrs_last, pos_last])
        data['hap2']['end'].append([chrs_last, pos_last])
        #for i in range(3, 8):
        #    pol = pollens[i-3]
        #    data[pol]['end'].append([chrs_last, pos_last])
 
    samples = ['hap1', 'hap2', 'Pollen41', 'Pollen42', 'Pollen43', 'Pollen44', 'Pollen45']
    fh = defaultdict(lambda : str())
    for sp in samples:
        fn = 'GT.{}.txt'.format(sp)
        ofile = open(fn, 'w')
        fh[sp] = ofile

    for sp in ['hap1', 'hap2']:
        for i in range(len(data[sp]['start'])):
            #print >> fh[sp], '{}\t{}\t{}\t{}'.format(data[sp]['start'][i][0], data[sp]['start'][i][1], data[sp]['end'][i][0], data[sp]['end'][i][1])
            #if sp == 'hap1':
            #    print data[sp]['start'][i][0]
            end_chrs = 'NA'
            end_pos  = 'NA'
            try:
                end_chrs = data[sp]['end'][i][0]
                end_pos  = data[sp]['end'][i][1]
            except:
                pass 
            #print >> fh[sp], '{}\t{}\t{}\t{}'.format(data[sp]['start'][i][0], data[sp]['start'][i][1], end_chrs, end_pos)
            #color = 'orange' if data[sp]['start'][i][2] == 1 else 'blue'
            color = 'gray'
            if data[sp]['start'][i][2] == 1:
                color = 'orange'
            elif data[sp]['start'][i][2] == 2:
                color = 'blue'
            print >> fh[sp], '{}\t{}\t{}\t{}\t{}\t+'.format(data[sp]['start'][i][0], data[sp]['start'][i][1], end_pos, data[sp]['start'][i][2], color)

    for sp in ['Pollen41', 'Pollen42', 'Pollen43', 'Pollen44', 'Pollen45']:
        for c in sorted(hap_chr_data[sp].keys()):
            for block in sorted(hap_chr_data[sp][c].keys(), key=int):
                start = np.min(hap_chr_data[sp][c][block])
                end   = np.max(hap_chr_data[sp][c][block])
                color = 'gray'
                if hap_chr_hap[sp][c][block] == 1:
                    color = 'orange'
                elif hap_chr_hap[sp][c][block] == 2:
                    color = 'blue'
                print >> fh[sp], '{}\t{}\t{}\t{}\t{}\t+'.format(c, start, end, hap_chr_hap[sp][c][block], color)

    for sp in samples:
        temp = fh[sp]
        temp.close()

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

    pollen_haplotype(args.input)

if __name__ == '__main__':
    main()

