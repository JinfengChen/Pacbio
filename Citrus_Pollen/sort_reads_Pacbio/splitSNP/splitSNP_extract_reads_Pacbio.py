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
1. Split list of fasta files
python splitSNP/splitSNP_extract_reads_Pacbio.py --input Pacbio_haplotype_reads --fasta_list input.fofn
2. Split fasta file
python splitSNP/splitSNP_extract_reads_Pacbio.py --input Pacbio_haplotype_reads --fasta citrus.pacbio.fasta

    '''
    print message


def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 60 --lines %s --interval 120 --task 1 --mem 15G --time 10:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)



def fasta_id(fastafile):
    fastaid = defaultdict(lambda : int())
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def fasta_id_list(fasta_list):
    fastaid = defaultdict(lambda : int())
    with open (fasta_list, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                fasta = line
                fastaid.update(fasta_id(fasta))
    return fastaid

def unique_reads_list(prefix, all_id):
    unique_alt_file = '%s.unique.alt.list' %(prefix)    
    unique_hom_ref_amb_file = '%s.unique.hom_ref_amb.list' %(prefix)
    if not os.path.exists(unique_alt_file):
        os.system('cat %s.alt.reads.list | sort | uniq > %s' %(prefix, unique_alt_file))
    if not os.path.exists(unique_hom_ref_amb_file):
    #    os.system('cat %s.hom.reads.list %s.ref.reads.list %s.amb.reads.list | sort | uniq > %s' %(prefix, prefix, prefix, unique_hom_ref_amb_file))
        ofile = open(unique_hom_ref_amb_file, 'w')
        alt_id = readtable(unique_alt_file)
        for id1 in sorted(all_id):
            if not alt_id.has_key(id1):
                print >> ofile, id1
        ofile.close()
    unique_hom_ref_amb_file = os.path.abspath(unique_hom_ref_amb_file)
    unique_alt_file         = os.path.abspath(unique_alt_file)
    return [unique_hom_ref_amb_file, unique_alt_file]

def split_file(script, fasta, hom_ref_amb, alt, prefix):
    hom_ref_amb_fasta = re.sub(r'.list', r'.fasta', hom_ref_amb)
    alt_fasta         = re.sub(r'.list', r'.fasta', alt)
    shell = '%s.fasta_split.sh' %(prefix)
    ofile = open(shell, 'w')
    print >> ofile, 'perl %s/getidseq.pl -l %s -f %s -o %s' %(script, hom_ref_amb, fasta, hom_ref_amb_fasta)
    print >> ofile, 'perl %s/getidseq.pl -l %s -f %s -o %s' %(script, alt, fasta, alt_fasta)
    ofile.close()
    runjob(shell, 1)

def split_file_list(script, fasta_list, hom_ref_amb, alt, prefix):
    dir_name_hom_ref_amb  = '%s_subreads_hom_ref_amb' %(prefix)
    dir_name_alt          = '%s_subreads_alt' %(prefix)
    if not os.path.exists(dir_name_hom_ref_amb):
        os.mkdir(dir_name_hom_ref_amb)
    if not os.path.exists(dir_name_alt):
        os.mkdir(dir_name_alt)
    shell = '%s.fasta_list_split.sh' %(prefix)
    ofile = open(shell, 'w')
    with open (fasta_list, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                fasta = line
                file_name = os.path.split(fasta)[1]
                hom_ref_amb_fasta = '%s/%s' %(os.path.abspath(dir_name_hom_ref_amb), file_name) 
                alt_fasta         = '%s/%s' %(os.path.abspath(dir_name_alt), file_name)
                print >> ofile, 'perl %s/getidseq.pl -l %s -f %s -o %s' %(script, hom_ref_amb, fasta, hom_ref_amb_fasta)
                print >> ofile, 'perl %s/getidseq.pl -l %s -f %s -o %s' %(script, alt, fasta, alt_fasta)
    ofile.close()
    runjob(shell, 2)

def readtable(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('--fasta', help='a fasta file to split')
    parser.add_argument('--fasta_list', help='a list of fasta file to split')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
 
    if not args.output:
        args.output = '%s_split_subreads' %(args.input)
    script = os.path.dirname(os.path.abspath(sys.argv[0]))

    #all sequence id
    all_fasta_id = defaultdict(lambda : int())
    if args.fasta:
        #all_fasta_id = fasta_id(args.fasta)
        all_fasta_id = readtable('citrus.pacbio.fasta.id') 
    if args.fasta_list:
        #all_fasta_id = fasta_id_list(args.fasta_list)  
        all_fasta_id = readtable('citrus.pacbio.fasta.id')    
    #unique id list
    hom_ref_amb, alt = unique_reads_list(args.input, all_fasta_id)
    #clear memory
    all_fasta_id.clear()
    #split fasta files
    if args.fasta:
        args.fasta = os.path.abspath(args.fasta)
        split_file(script, args.fasta, hom_ref_amb, alt, args.input)
    if args.fasta_list:
        args.fasta_list = os.path.abspath(args.fasta_list)
        split_file_list(script, args.fasta_list, hom_ref_amb, alt, args.input)

if __name__ == '__main__':
    main()

