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

def get_gene_mrna_id(gff):
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
            print mRNA.id

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

    gene_info = get_gene_mrna_id(args.gff)    

if __name__ == '__main__':
    main()


