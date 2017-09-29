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
python Prepare_10x.py --vcf phased_variants.vcf.gz --haplotype Haplotype.FCM.txt

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

'''
chr1:557344	C	A	A	A	C	NA	A	C	A
chr1:557377	T	A	A	A	T	NA	A	T	A
'''
def read_haplotype(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = [unit[8], unit[9]]
    return data

'''
chr1    557344  .       C       A       2349.72 PASS    MUMAP_REF=59.86;MUMAP_ALT=60.0;AO=44;RO=50;MMD=1.325288775;RESCUED=0;NOT_RESCUED=94     GT:BX:PS:PQ:JQ  0|1:
'''
def summary_haplotype_block(vcf, haplotype, outfile_up, outfile_down):
    phased_block = defaultdict(lambda : list())
    phased_block_com  = defaultdict(lambda : list())
    phased_block_rank = defaultdict(lambda : int())
    count = 0
    import VCF
    for v in VCF.lines(vcf):
        block_id = re.split(r':', v['FCM'][1])[1]
        genotype = re.split(r':', v['FCM'][0])[0]
        #print genotype
        #print '{}\t{}\t{}\t{}'.format(v['CHROM'], v['POS'], block_id, v['FCM'][0])
        block_idx = '{}_{}'.format(v['CHROM'], block_id)
        snp_idx   = '{}:{}'.format(v['CHROM'], v['POS'])
        bases     = [v['REF'], v['ALT']]
        hap1_10x  = bases[int(genotype[0])]
        haplotype_flag = -1
        if haplotype.has_key(snp_idx):
            if hap1_10x == haplotype[snp_idx][0]:
                haplotype_flag = 0
            elif hap1_10x == haplotype[snp_idx][1]:
                haplotype_flag = 1
        if not phased_block.has_key(block_idx):
            count += 1
        phased_block[block_idx].append(int(v['POS']))
        phased_block_com[block_idx].append(haplotype_flag)
        phased_block_rank[block_idx] = count

    for blc in phased_block.keys():
        chrs, blc_id = re.split(r'_', blc)
        snps  = len(phased_block[blc])  
        start = np.min(phased_block[blc])
        end   = np.max(phased_block[blc])
        length= int(end) - int(start) + 1
        hap1_n = len([i for i in phased_block_com[blc] if i == 0])
        hap2_n = len([i for i in phased_block_com[blc] if i == 1])
        hap0_n = len([i for i in phased_block_com[blc] if i == -1])
        print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(blc, snps, length, chrs, start, end, hap1_n, hap2_n, hap0_n)
        ratio = 0.0
        #try:
        #    ratio = np.max(float(hap1_n)/(float(hap1_n)+float(hap2_n)), float(hap2_n)/(float(hap1_n)+float(hap2_n)))
        #except:
        #    continue
        min_snp = 1
        if hap1_n >= min_snp and hap1_n == snps - hap0_n:
            ratio = 1
        elif hap2_n >= min_snp and hap2_n == snps - hap0_n:
            ratio = 1
        elif hap1_n >= min_snp and hap1_n > hap2_n:
            ratio = float(hap1_n)/(float(hap1_n)+float(hap2_n))
        elif hap2_n >= min_snp and hap2_n > hap1_n:
            ratio = float(hap2_n)/(float(hap1_n)+float(hap2_n)) 

        print ratio
        if ratio > 0.95: 
            color = 'gray'
            hap   = 3
            print hap1_n, hap2_n
            if hap1_n > hap2_n:
                color = 'orange'
                hap   = 1
            else:
                color = 'blue' 
                hap   = 2
            print color
            if phased_block_rank[blc]%2 == 1:
                print >> outfile_up, '{}\t{}\t{}\t{}\t{}\t+'.format(chrs, start, end, hap, color)
            else:
                print >> outfile_down, '{}\t{}\t{}\t{}\t{}\t+'.format(chrs, start, end, hap, color)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf')
    parser.add_argument('--haplotype')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    try:
        len(args.vcf) > 0 and len(args.haplotype) > 0
    except:
        usage()
        sys.exit(2)


    pollen_haplotype = read_haplotype(args.haplotype)

    outfile_up   = open('GT.10x_up.unsort.txt', 'w')
    outfile_down = open('GT.10x_down.unsort.txt', 'w')
    summary_haplotype_block(args.vcf, pollen_haplotype, outfile_up, outfile_down)
    outfile_up.close()
    outfile_down.close()
   
    os.system('sort -k1,1 -k2,2n GT.10x_up.unsort.txt > GT.10x_up.txt')       
    os.system('sort -k1,1 -k2,2n GT.10x_down.unsort.txt > GT.10x_down.txt') 

if __name__ == '__main__':
    main()

