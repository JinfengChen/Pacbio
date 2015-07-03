#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
#import numpy as np
import re
import os
import argparse
from Bio import SeqIO
from pbcore.io import AlignmentSet
from pbcore.io import ReferenceSet
from GenomicConsensus import reference
import argparse, atexit, cProfile, gc, glob, h5py, logging, multiprocessing
import os, pstats, random, shutil, tempfile, time, threading, Queue, traceback


def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


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
    parser.add_argument('-r', '--ref')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    with AlignmentSet(args.input) as peekCmpH5:
        print "test"
        print peekCmpH5
        logging.info("Peeking at CmpH5 file %s" % (args.input))
        logging.info("Input CmpH5 data: numAlnHits=%d" % len(peekCmpH5))
        cmpContigNames = set(peekCmpH5.refNames)
        print cmpContigNames
        reference.loadFromFile(args.ref, peekCmpH5)
        
    f = ReferenceSet(args.ref)
    f.assertIndexed()
    for fastaRecord in f.contigs:
        refName = fastaRecord.id
        print refName
        #if refName in cmpContigNames:



if __name__ == '__main__':
    main()

