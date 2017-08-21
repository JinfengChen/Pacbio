#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO

import gffutils

def usage():
    test="name"
    message='''
usage: long_intron_gene.py [-h] [--gff GFF] [-o OUTPUT] [-v]

optional arguments:
  -h, --help            show this help message and exit
  --gff GFF             gene gff
  -o OUTPUT, --output OUTPUT
  -v

    '''
    print message

def get_gene_exon_num(gff):
    if not os.path.exists('{}.db'.format(gff)):
        gffutils.create_db(gff, dbfn='{}.db'.format(gff), merge_strategy="create_unique")
    gff_db = gffutils.FeatureDB('{}.db'.format(gff), keep_order=True)
    gene_exon_num = defaultdict(lambda : int())
    for gene in gff_db.features_of_type('gene', order_by='start'):
        mRNA_longest = ''
        mRNA_len     = 0
        utr_flag     = 0
        utr_mrna_t   = 0
        utr_long_t   = 0
        longest_intron = 0
        for mRNA in gff_db.children(gene, featuretype="mRNA"):
            cds_len = 0
            last_exon_end = 0
            for CDS in gff_db.children(mRNA, featuretype="CDS", order_by='start'):
                cds_len += CDS[4] - CDS[3] + 1
                if last_exon_end == 0:
                   last_exon_end = CDS[4]
                   continue
                else:
                   current_intron = abs(CDS[3] - last_exon_end)
                   #if gene.id == "evm.TU.chr1.101":
                   #    print '{} {} {} {} {} {} {}'.format(gene, mRNA, CDS, CDS[3], last_exon_end, current_intron, longest_intron)
  
                   if current_intron > longest_intron:
                       longest_intron = current_intron
                   last_exon_end = CDS[4]
            if cds_len > mRNA_len:
                mRNA_len     = cds_len
                mRNA_longest = mRNA.id 
            utr_mrna_t += 1
            for UTR in gff_db.children(mRNA, featuretype=['three_prime_UTR', 'five_prime_UTR']):
                utr_len = UTR[4] - UTR[3] + 1
                if utr_len > 1500:
                    utr_long_t += 1
        if utr_long_t == utr_mrna_t:
            utr_flag = 1
        exon_num = len(list(gff_db.children(mRNA_longest, featuretype='CDS', level=1)))
        exon_len = 0
        for exon in list(gff_db.children(mRNA_longest, featuretype='CDS', level=1)):
            exon_len += exon[4] - exon[3] + 1
        gene_exon_num[gene.id] = [exon_num, exon_len, utr_flag, longest_intron, gene]
        #print('{}\t{}'.format(gene.id, exon_num))
    return gene_exon_num


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', help="gene gff")
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)

    gene_info = get_gene_exon_num(args.gff)    

    for gene in sorted(gene_info.keys()):
        if gene_info[gene][3] > 10000:
            print '{}: {} {}'.format(gene, gene_info[gene][3], gene_info[gene][4])

if __name__ == '__main__':
    main()


