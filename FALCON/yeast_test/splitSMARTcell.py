#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

#m131225_223149_42137_c100577112550000001823096901191465_s1_p0/108990/0_5501
def fasta_id(fastafile):
    r = re.compile(r'(.*?)\/')
    file_d = []
    for record in SeqIO.parse(fastafile, "fasta"):
        if r.search(str(record.id)):
            ofile = open('%s.subread.fa' %(r.search(str(record.id)).groups(0)[0]), 'a')
            file_d.append(ofile)
            print >> ofile, '>%s\n%s' %(str(record.id), str(record.seq))
            ofile.close()
    #for f in file_d:
    #    f.close()


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
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    fasta_id(args.input)

if __name__ == '__main__':
    main()

