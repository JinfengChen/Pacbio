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
import subprocess


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

def update_id(gff, ids, rewrite):
    if not os.path.exists('{}.db'.format(gff)):
        gffutils.create_db(gff, dbfn='{}.db'.format(gff), merge_strategy="create_unique")
    gff_db = gffutils.FeatureDB('{}.db'.format(gff), keep_order=True)
    
    #gff_out_file_temp = '{}.locus_tag_id.temp.gff'.format(os.path.splitext(gff)[0])
    gff_out_file      = '{}.locus_tag_id.gff'.format(os.path.splitext(gff)[0])
    
    if os.path.exists(gff_out_file) and rewrite == 0:
        return gff_out_file
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
            mRNA_rec_id = '{}.{}'.format(ids[gene.id], mRNA_count)
            #mRNA_attrs_raw = dict(mRNA_rec.attributes)
            #for k in mRNA_attrs.keys():
            #    if not k in ['ID', 'Parent']:
            #        del mRNA_attrs[k]
            #mRNA_attrs['ID'] = [mRNA_rec_id]
            #mRNA_attrs['Parent'] = [ids[gene.id]]
            mRNA_attrs = {'ID':[mRNA_rec_id], 'Parent':[ids[gene.id]], 'Name':[attrs['Name'][0]]}
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
            write_mRNA_children(gff_out, gff_db, mRNA_rec, mRNA_rec_id)
    #clean_cmd = "sed 's/5_prime_partial=true;//g' {} | sed 's/3_prime_partial=true;//g' > {}".format(gff_out_file_temp, gff_out_file)
    #print clean_cmd
    #os.system(clean_cmd)
    #failure = subprocess.call(clean_cmd, shell=True)
    return gff_out_file

def write_mRNA_children(gff_out, db, mRNA_rec, mRNA_rec_id):
    print 'writing mRNA childern: {}'.format(mRNA_rec_id)
    mRNA_children = db.children(mRNA_rec, order_by='start')
    features = ['exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'] 
    for child_rec in mRNA_children:
        print child_rec
        counter = defaultdict(lambda : int())
        if child_rec.featuretype in features:
            ft = child_rec.featuretype
            counter[ft] += 1
            #child_rec.id = '{}:exon{}'.format(mRNA_rec.id, exon_count)
            attrs = dict(child_rec.attributes)
            attrs['ID'] = ['{}:{}{}'.format(mRNA_rec_id, ft, counter[ft])]
            attrs['Parent'] = [mRNA_rec_id]
            child_rec = gffutils.Feature(
                seqid=child_rec.chrom,
                source='Fairchild_v1.0',
                featuretype='{}'.format(ft),
                start=child_rec.start,
                end=child_rec.end,
                strand=child_rec.strand,
                attributes=attrs)
            gff_out.write_rec(child_rec)

def write_pep_from_gff(gff, genome, tools):
    prefix = os.path.splitext(gff)[0]
    cmds = []
    cmds.append('{} {} {} > {}.cds.fa'.format(tools['getgene'], gff, genome, prefix))
    cmds.append('{} {}.cds.fa > {}.pep.fa'.format(tools['cds2aa'], prefix, prefix))
    if os.path.exists('{}.pep.fa'.format(prefix)):
        return
    for cmd in cmds:
        print(cmd)
        os.system(cmd)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff')
    parser.add_argument('--id')
    parser.add_argument('--genome')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0 and len(args.id) > 0
    except:
        usage()
        sys.exit(2)
   
    tools = defaultdict(lambda : str())
    script_path = '{}/scripts'.format(os.path.dirname(os.path.realpath(sys.argv[0]))) 
    tools["solar"]     = 'perl ~/BigData/00.RD/Annotation/HEG4/protein-map-genome/bin/solar/solar.pl'
    tools["bestalign"] = 'perl {}/bestAlign.pl'.format(script_path)
    tools["fastadeal"] = 'perl {}/fastaDeal.pl'.format(script_path)
    tools["getidseq="] = 'perl {}/getidseq.pl'.format(script_path)
    tools["getidgff"]  = 'perl {}/getidgff.pl'.format(script_path)
    tools["getgene"]  = 'perl {}/getGene.pl'.format(script_path)   
    tools["cds2aa"]  = 'perl {}/cds2aa.pl'.format(script_path)
    tools["blast2blastm8"]  = 'perl {}/blast2blastm8.pl'.format(script_path)


    ids = defaultdict(lambda : str())
    for line in np.loadtxt(args.id, dtype=str):
        ids[line[0]] = line[1]
    new_gff = update_id(args.gff, ids, 0)
    write_pep_from_gff(new_gff, args.genome, tools) 

if __name__ == '__main__':
    main()

