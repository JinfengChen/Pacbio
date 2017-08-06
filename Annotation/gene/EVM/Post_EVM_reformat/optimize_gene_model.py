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
usage: optimize_gene_model.py [-h] [--evm_gff EVM GFF] [--maker_gff MAKER GFF]
                              [--pasa_as_gff PASA ALTERNATIVE SPLICING GFF]
                              [--repeat_gff TE GFF]
                              [--TE_blastp EVM PEP TO TE PROTEIN BLASTP RESULTS]
                              [--US_blastp EVM PEP TO UNIPROT_SPROT PROTEIN BLASTP RESULTS]
                              [-o OUTPUT] [-v]

optional arguments:
  -h, --help            show this help message and exit
  --evm_gff EVM GFF
  --maker_gff MAKER GFF
  --pasa_as_gff PASA ALTERNATIVE SPLICING GFF
  --repeat_gff TE GFF
  --TE_blastp EVM PEP TO TE PROTEIN BLASTP RESULTS
  --US_blastp EVM PEP TO UNIPROT_SPROT PROTEIN BLASTP RESULTS
  -o OUTPUT, --output OUTPUT
  -v



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

def filter_repeat_gene(gene_gff, repeat_gff, gene_swissprot_id, gene_exon_num):
    from pybedtools import BedTool
    repeats = BedTool(repeat_gff) 
    genes   = BedTool(gene_gff)
    genes_in_repeat = defaultdict(lambda : int())
    for g in genes.intersect(repeats):
        if g[2] == 'gene':
                temp  = defaultdict(str)
                attrs = re.split(r';', g[8])
                for attr in attrs:
                    if not attr == '':
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                gene_id   = temp['ID']
                genes_in_repeat[gene_id] = 1
    genes_low_quality = defaultdict(lambda : int())
    count_temp = 0
    for g in sorted(genes_in_repeat.keys()): 
        if not gene_swissprot_id.has_key(g):
            count_temp += 1
            if gene_exon_num.has_key(g):
                if gene_exon_num[g][0] <= 2 and gene_exon_num[g][1] <= 600:
                    genes_low_quality[g] = 1
    print 'number of gene in repeat {}'.format(len(genes_in_repeat.keys()))
    print 'number of gene in repeat without swissprot hit {}'.format(count_temp)
    print 'number of gene i low quality {}'.format(len(genes_low_quality.keys()))
    return genes_low_quality
 
def write_pep_from_gff(gff, genome, tools):
    prefix = os.path.splitext(gff)[0]
    cmds = []
    cmds.append('{} {} {} > {}.cds.fa'.format(tools['getgene'], gff, genome, prefix))
    cmds.append('{} {}.cds.fa > {}.pep.fa'.format(tools['cds2aa'], prefix, prefix))
    for cmd in cmds:
        print(cmd)
        os.system(cmd)
 


def write_gene_gff(gff, gene_repeat_id, prefix, subtitle, reverse):
    if not os.path.exists('{}.db'.format(gff)):
        gffutils.create_db(gff, dbfn='{}.db'.format(gff), merge_strategy="create_unique")
    gff_db = gffutils.FeatureDB('{}.db'.format(gff), keep_order=True)
    gff_noTE_file = '{}.{}.gff'.format(prefix, subtitle)
    if os.path.exists(gff_noTE_file):
        return gff_noTE_file

    gff_noTE = gffutils.gffwriter.GFFWriter(gff_noTE_file)
    for gene in gff_db.features_of_type('gene', order_by='start'):
        if reverse:
            if not gene_repeat_id.has_key(gene.id):
                gff_noTE.write_gene_recs(gff_db, gene.id)
        else:
            if gene_repeat_id.has_key(gene.id):
                gff_noTE.write_gene_recs(gff_db, gene.id)
    return gff_noTE_file

def get_gene_exon_num(gff):
    if not os.path.exists('{}.db'.format(gff)):
        gffutils.create_db(gff, dbfn='{}.db'.format(gff), merge_strategy="create_unique")
    gff_db = gffutils.FeatureDB('{}.db'.format(gff), keep_order=True)
    gene_exon_num = defaultdict(lambda : int())
    for gene in gff_db.features_of_type('gene', order_by='start'):
        exon_num = len(list(gff_db.children(gene, featuretype='CDS', level=2)))
        exon_len = 0
        for exon in list(gff_db.children(gene, featuretype='CDS', level=2)):
            exon_len += exon[4] - exon[3] + 1
        gene_exon_num[gene.id] = [exon_num, exon_len]
        #print('{}\t{}'.format(gene.id, exon_num))
    return gene_exon_num

def convert_mrna2gene(mrna, gff):
    if not os.path.exists('{}.db'.format(gff)):
        gffutils.create_db(gff, dbfn='{}.db'.format(gff), merge_strategy="create_unique")
    gff_db = gffutils.FeatureDB('{}.db'.format(gff), keep_order=True)
    gene_id = defaultdict(lambda : int())
    for mrna_id in sorted(mrna):
        gene = [g.id for g in gff_db.parents(mrna_id, featuretype="gene", level=1)]
        gene_id[gene[0]] = 1
        #print '{}: {}'.format(mrna_id, gene[0])
    return gene_id

def blast_best_hit_list_us(blast, tools, cutoff):
    if not os.path.exists(blast):
        print '{} not exists'.format(blast)
        sys.exit(2)
    table_out = '{}.m8'.format(os.path.splitext(blast)[0])
    cmds = []
    if not os.path.exists(table_out):
        cmds.append('{} {} > {}'.format(tools["blast2blastm8"], blast, table_out))
    if not os.path.exists(solar_out):
        cmds.append('{} -d -1 {} > {}'.format(tools["solar"], table_out, solar_out)) 
    if len(cmds) > 0:
        for cmd in cmds:
            print cmd
            os.system(cmd)
    

def blast_best_hit_list(blast, tools, cutoff):
    if not os.path.exists(blast):
        print '{} not exists'.format(blast)
        sys.exit(2)

    solar_out = '{}.solar'.format(os.path.splitext(blast)[0])
    best_id   = '{}.solar.id'.format(os.path.splitext(blast)[0])
    cmds = []
    if not os.path.exists(solar_out):
        cmds.append('{} -d -1 {} > {}'.format(tools["solar"], blast, solar_out))
    if not os.path.exists(best_id):
        cmds.append('{} {} -cutoff {} | cut -f1 > {}'.format(tools["bestalign"], solar_out, cutoff, best_id))
    if len(cmds) > 0:
        for cmd in cmds:
            print cmd
            os.system(cmd)
    return best_id

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome', help="genome assembly")
    parser.add_argument('--evm_gff', help="EVM gff")
    parser.add_argument('--maker_gff', help="Maker gff")
    parser.add_argument('--pasa_as_gff', help="PASA alternative splicing gff")
    parser.add_argument('--repeat_gff', help="TE gff")
    parser.add_argument('--TE_blastp', help="evm pep to TE protein blastp results")
    parser.add_argument('--US_blastp', help="evm pep to uniprot_sprot protein blastp results")
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.evm_gff) > 0
        len(args.maker_gff) > 0
        len(args.pasa_as_gff) > 0
        len(args.repeat_gff) > 0
        len(args.TE_blastp) > 0
        len(args.US_blastp) > 0
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

 
    #filter by Tpase
    mrna_repeat_id_file = blast_best_hit_list(args.TE_blastp, tools, 0.7)
    mrna_repeat_id      = list(np.loadtxt(mrna_repeat_id_file, dtype=str)) 
    gene_repeat_id = convert_mrna2gene(mrna_repeat_id, args.evm_gff) 
    evm_noTE_gff   = write_gene_gff(args.evm_gff, gene_repeat_id, args.output, 'noTE', 1)
    write_pep_from_gff(evm_noTE_gff, args.genome, tools)

    #filter by repeatmasker and swisspro
    mrna_swissprot_id_file = blast_best_hit_list(args.US_blastp, tools, 0.3) 
    mrna_swissprot_id      = list(np.loadtxt(mrna_swissprot_id_file, dtype=str))
    gene_swissprot_id      = convert_mrna2gene(mrna_swissprot_id, args.evm_gff)
    gene_exon_num          = get_gene_exon_num(args.evm_gff) 
    gene_low_quality_id    = filter_repeat_gene(args.evm_gff, args.repeat_gff, gene_swissprot_id, gene_exon_num)
    evm_noTE_highqual_gff  = write_gene_gff(evm_noTE_gff, gene_low_quality_id, args.output, 'noTE_highqual', 1)
    write_pep_from_gff(evm_noTE_highqual_gff, args.genome, tools)

if __name__ == '__main__':
    main()

