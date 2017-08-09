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
import gffutils


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

def update_id(gff, ids):
    if not os.path.exists('{}.db'.format(gff)):
        gffutils.create_db(gff, dbfn='{}.db'.format(gff), merge_strategy="create_unique")
    gff_db = gffutils.FeatureDB('{}.db'.format(gff), keep_order=True)
    
    gff_out_file = '{}.locus_tag_id.gff'.format(os.path.splitext(gff)[0])
    gff_out = gffutils.gffwriter.GFFWriter(gff_out_file, with_header=False) 
    for gene in gff_db.features_of_type('gene', order_by=['seqid', 'start']):
        #write gene
        #gene_rec = gff_db[gene]
        #print '{} : {} : {}'.format(gene_rec.id, gene.id, ids[gene.id])
        attrs = dict(gene.attributes)
        attrs['ID'] = [ids[gene.id]]
        gene_rec = gffutils.Feature(
            seqid=gene.chrom,
            source='Fairchild_v1.0',
            featuretype='gene',
            start=gene.start,
            end=gene.end,
            strand=gene.strand,
            attributes=attrs)
        gff_out.write_rec(gene_rec)

        #write mRNA
        mRNA_lens = {}
        c = list(gff_db.children(gene, featuretype="mRNA"))
        for mRNA in gff_db.children(gene, featuretype="mRNA"):
            mRNA_lens[mRNA.id] = \
                sum(len(exon) for exon in gff_db.children(mRNA,
                                                      featuretype="exon"))
        sorted_mRNAs = \
            sorted(mRNA_lens.items(), key=lambda x: x[1], reverse=True)
        mRNA_count = 0
        for curr_mRNA in sorted_mRNAs:
            mRNA_count += 1
            mRNA_id = curr_mRNA[0]
            mRNA_rec = gff_db[mRNA_id]
            #mRNA_rec.id = '{}.{}'.format(ids[gene], mRNA_count)
            mRNA_attrs = dict(mRNA_rec.attributes)
            mRNA_attrs['ID'] = ['{}.{}'.format(ids[gene.id], mRNA_count)]
            mRNA_attrs['Parent'] = ids[gene.id]
            mRNA_rec_new = gffutils.Feature(
                seqid=mRNA_rec.chrom,
                source='Fairchild_v1.0',
                featuretype='mRNA',
                start=mRNA_rec.start,
                end=mRNA_rec.end,
                strand=mRNA_rec.strand,
                attributes=mRNA_attrs)
            gff_out.write_rec(mRNA_rec_new)
            #write exon/cds/utr
            write_mRNA_children(gff_out, gff_db, mRNA_rec, mRNA_rec_new) 


def write_mRNA_children(gff_out, db, mRNA_rec_old, mRNA_rec):
    mRNA_children = db.children(mRNA_rec_old.id, order_by='start')
    features = ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'] 
    for child_rec in mRNA_children:
        counter = defaultdict(lambda : int())
        for ft in features:
            counter[ft] += 1
            #child_rec.id = '{}:exon{}'.format(mRNA_rec.id, exon_count)
            attrs = dict(child_rec.attributes)
            attrs['ID'] = ['{}:{}{}'.format(mRNA_rec.id, ft, counter[ft])]
            child_rec = gffutils.Feature(
                seqid=child_rec.chrom,
                source='Fairchild_v1.0',
                featuretype='{}'.format(ft),
                start=child_rec.start,
                end=child_rec.end,
                strand=child_rec.strand,
                attributes=attrs)
            gff_out.write_rec(child_rec)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff')
    parser.add_argument('--id')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0 and len(args.id) > 0
    except:
        usage()
        sys.exit(2)

    ids = defaultdict(lambda : str())
    for line in np.loadtxt(args.id, dtype=str):
        ids[line[0]] = line[1]
    update_id(args.gff, ids)

if __name__ == '__main__':
    main()

