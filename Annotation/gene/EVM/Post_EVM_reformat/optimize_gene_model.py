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
    mrna_info = defaultdict(lambda : list())
    for gene in sorted(gene_exon_num.keys()):
        mrna   = gene_exon_num[gene][4]
        cdsn   = gene_exon_num[gene][0]
        #print 'longest mRNA: {} {} {}'.format(mrna, gene, cdsn)
        mrna_info[mrna] = [gene, cdsn]

    from pybedtools import BedTool
    repeats = BedTool(repeat_gff) 
    genes   = BedTool(gene_gff)
    genes_in_repeat = defaultdict(lambda : int())
    genes_in_repeat_nonprotein = defaultdict(lambda : int())
    mrna_repeat_dict = defaultdict(lambda : defaultdict(lambda : int()))
    for g in genes.intersect(repeats):
        #print g
        if g[2] == 'CDS':
                temp  = defaultdict(str)
                attrs = re.split(r';', g[8])
                for attr in attrs:
                    if not attr == '':
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                cds_id   = '{}_{}_{}'.format(temp['ID'], g[3], g[4])
                mrna_id  = temp['Parent']
                #print 'CDS: {}, {}'.format(mrna_id, cds_id)
                mrna_repeat_dict[mrna_id][cds_id] = 1
    for m in sorted(mrna_repeat_dict.keys()):
        cdsn = len(mrna_repeat_dict[m].keys())
        #print 'mRNA: {}, {}'.format(m, cdsn)
        if mrna_info.has_key(m):
            if float(cdsn)/mrna_info[m][1] > 0.7:
                gene_id = mrna_info[m][0]
                genes_in_repeat[gene_id] = 1
    genes_low_quality = defaultdict(lambda : int())
    count_temp = 0
    for g in sorted(genes_in_repeat.keys()): 
        if not gene_swissprot_id.has_key(g):
            count_temp += 1
            genes_in_repeat_nonprotein[g] = 1
            #if gene_exon_num.has_key(g):
                #if gene_exon_num[g][0] <= 2 and gene_exon_num[g][1] <= 600:
                #    genes_low_quality[g] = 1
            genes_low_quality[g] = 1    
    print 'number of gene in repeat {}'.format(len(genes_in_repeat.keys()))
    print 'number of gene in repeat without swissprot hit {}'.format(count_temp)
    print 'number of gene i low quality {}'.format(len(genes_low_quality.keys()))
    return genes_low_quality, genes_in_repeat
 
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
 
def replace_gene_fusion(gene_info, maker_gff_db):
    print 'in replace function: {}, {}, {}'.format(gene_info.seqid, gene_info.start, gene_info.end)
    predictor_num = defaultdict(lambda : list())
    for overlap_gene in maker_gff_db.region(region=(gene_info.seqid, gene_info.start, gene_info.end), featuretype=['gene']):
        #print overlap_gene
        if overlap_gene[1] == 'maker':
            predictor_num['maker'].append(overlap_gene.id)
        elif overlap_gene[1] == 'fgenesh_masked':
            predictor_num['fgenesh'].append(overlap_gene.id)
        elif overlap_gene[1] == 'augustus_masked':
            predictor_num['augustus'].append(overlap_gene.id)
        elif overlap_gene[1] == 'snap_masked':
            predictor_num['snap'].append(overlap_gene.id)
        elif overlap_gene[1] == 'genemark':
            predictor_num['genemark'].append(overlap_gene.id)
    fusion_count = 0
    for p in predictor_num.keys():
        if len(predictor_num[p]) == 2:
            fusion_count += 1
    #print 'fusion: {}'.format(fusion_count)
    if fusion_count >= 3:
        print 'fusion: {}'.format(fusion_count)
        print 'gene info: {}, {}, {}'.format(gene_info.seqid, gene_info.start, gene_info.end)
        if len(predictor_num['maker']) == 2:
            return predictor_num['maker']
        elif len(predictor_num['fgenesh']) == 2:
            return predictor_num['fgenesh']
        elif len(predictor_num['augustus']) == 2:            
            return predictor_num['augustus']
        elif len(predictor_num['genemark']) == 2:
            return predictor_num['genemark']
        else:
            return ['']
    else: 
        return ['']

def read_gene_fusion_list(gene_fusion_list):
    gene_fusion_dict = defaultdict(lambda : list())
    for line in np.loadtxt(gene_fusion_list, dtype=str): 
        evm_genes   = re.split(r';', line[0])
        maker_genes = re.split(r';', line[1])
        for i in range(len(evm_genes)):
            if i == 0:
                gene_fusion_dict[evm_genes[i]] = maker_genes
            else:
                gene_fusion_dict[evm_genes[i]] = ['']
    return gene_fusion_dict

def select_gene_model(source_gene_id, source_gff, pasa_as_gff, maker_gff, prefix, subtitle, gene_fusion_list, gene_in_repeat_id):
    if not os.path.exists('{}.db'.format(source_gff)):
        gffutils.create_db(source_gff, dbfn='{}.db'.format(source_gff), merge_strategy="create_unique")
    source_gff_db = gffutils.FeatureDB('{}.db'.format(source_gff), keep_order=True)

    if not os.path.exists('{}.db'.format(pasa_as_gff)):
        gffutils.create_db(pasa_as_gff, dbfn='{}.db'.format(pasa_as_gff), merge_strategy="create_unique")
    pasa_as_gff_db = gffutils.FeatureDB('{}.db'.format(pasa_as_gff), keep_order=True)

    if not os.path.exists('{}.db'.format(maker_gff)):
        gffutils.create_db(maker_gff, dbfn='{}.db'.format(maker_gff), merge_strategy="create_unique")
    maker_gff_db = gffutils.FeatureDB('{}.db'.format(maker_gff), keep_order=True)

    source_gene_exons  = get_gene_exon_num(source_gff)
    pasa_as_gene_exons = get_gene_exon_num(pasa_as_gff)   

    gff_out_file = '{}.{}.gff'.format(prefix, subtitle)
    if os.path.exists(gff_out_file):
        return gff_out_file

    gene_fusion_dict = read_gene_fusion_list(gene_fusion_list)

    intron_cutoff  = 20000
    count_total    = 0
    count_pasa     = 0
    count_fusion   = 0
    count_utr_long = 0
    count_intron_long  = 0
    count_intron_long2 = 0
    count_no_pasa  = 0
    count_exon_num = 0
    count_exon_len = 0
    gff_out = gffutils.gffwriter.GFFWriter(gff_out_file)    
    for gene in sorted(source_gene_id.keys()):
        print 'process gene model: {}'.format(gene)
        count_total += 1

        '''gene_fusion_model = replace_gene_fusion(source_gff_db[gene], maker_gff_db)
        print 'fusion model: {}'.format(gene_fusion_model)
        if len(gene_fusion_model) > 1:
            count_fusion += 1
            print '{}: {}'.format(gene, source_gff_db[gene])
            print 'replaced: {}'.format(gene_fusion_model)
            for temp_gene in sorted(gene_fusion_model):
                gff_out.write_gene_recs(maker_gff_db, temp_gene)
            continue
        '''
        if gene_fusion_dict.has_key(gene):
            count_fusion += 1
            print '{}: {}'.format(gene, source_gff_db[gene])
            print 'replaced: {}'.format(gene_fusion_dict[gene])
            for temp_gene in gene_fusion_dict[gene]:
                if not temp_gene == '' and not temp_gene == 'NA':
                    gff_out.write_gene_recs(maker_gff_db, temp_gene)
            continue

        if not pasa_as_gene_exons.has_key(gene):
            count_no_pasa += 1
            if source_gene_exons[gene][3] > intron_cutoff and gene_in_repeat_id.has_key(gene):
                count_intron_long += 1
                print 'long intron gene out: {}'.format(gene)
                continue
            else:
                if source_gene_exons[gene][3] > intron_cutoff:
                    print 'long intron gene still in: {}'.format(gene)
                    count_intron_long2 += 1
                gff_out.write_gene_recs(source_gff_db, gene)
                continue

        if source_gene_exons[gene][0] == pasa_as_gene_exons[gene][0]:
            #if pasa_as_gene_exons[gene][2] == 1:
            #    count_utr_long += 1
            #    gff_out.write_gene_recs(source_gff_db, gene)
            #    continue
            if pasa_as_gene_exons[gene][3] > intron_cutoff and gene_in_repeat_id.has_key(gene):
                count_intron_long += 1
                #gff_out.write_gene_recs(source_gff_db, gene)
                print 'long intron gene out: {}'.format(gene)
                continue
            elif pasa_as_gene_exons[gene][3] > intron_cutoff:
                print 'long intron gene still in: {}'.format(gene)
                count_intron_long2 += 1

            if pasa_as_gene_exons[gene][2] == 1:
                count_utr_long += 1
                gff_out.write_gene_recs(source_gff_db, gene)
                continue

            if pasa_as_gene_exons[gene][1] >= 0.9 * source_gene_exons[gene][1] and pasa_as_gene_exons[gene][1] <= 1.1 * source_gene_exons[gene][1]:
                count_pasa += 1
                gff_out.write_gene_recs(pasa_as_gff_db, gene)
            else:
                count_exon_len += 1
                gff_out.write_gene_recs(source_gff_db, gene)
        else:
            if source_gene_exons[gene][3] > intron_cutoff and gene_in_repeat_id.has_key(gene):
                count_intron_long += 1
                print 'long intron gene out: {}'.format(gene)
                continue
            else:
                if source_gene_exons[gene][3] > intron_cutoff:
                    print 'long intron gene still in: {}'.format(gene)
                    count_intron_long2 += 1
                #gff_out.write_gene_recs(source_gff_db, gene)
                #count_exon_num += 1
                #continue 
            count_exon_num += 1
            gff_out.write_gene_recs(source_gff_db, gene)

    print 'total gene model: {}'.format(count_total)
    print 'genes without pasa models: {}'.format(count_no_pasa)
    print 'gene fusion: {}'.format(count_fusion)
    print 'genes utr too long (>1.5 kb): {}'.format(count_utr_long)
    print 'genes intron too long out (>10 kb): {}'.format(count_intron_long)
    print 'genes intron too long in (>10 kb): {}'.format(count_intron_long2)
    print 'genes cds number different: {}'.format(count_exon_num)
    print 'genes cds length different: {}'.format(count_exon_len)
    print 'genes write as pasa model: {}'.format(count_pasa)
    return gff_out_file 

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

def get_gene_id(gff):
    if not os.path.exists('{}.db'.format(gff)):
        gffutils.create_db(gff, dbfn='{}.db'.format(gff), merge_strategy="create_unique")
    gff_db = gffutils.FeatureDB('{}.db'.format(gff), keep_order=True)
    gene_id = defaultdict(lambda : int())
    for gene in gff_db.features_of_type('gene', order_by='start'):
        gene_id[gene.id] = 1
    return gene_id 

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
        gene_exon_num[gene.id] = [exon_num, exon_len, utr_flag, longest_intron, mRNA_longest]
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

'''
def blast_best_hit_list_us(blast, tools, cutoff):
    if not os.path.exists(blast):
        print '{} not exists'.format(blast)
        sys.exit(2)
    cmds = []
    table_out = '{}.m8'.format(os.path.splitext(blast)[0])
    if blast.endswith('.m8'):
        table_out = blast
    else: 
        if not os.path.exists(table_out):
            cmds.append('{} {} > {}'.format(tools["blast2blastm8"], blast, table_out))
    if not os.path.exists(solar_out):
        cmds.append('{} -d -1 {} > {}'.format(tools["solar"], table_out, solar_out)) 
    if len(cmds) > 0:
        for cmd in cmds:
            print cmd
            os.system(cmd)
'''    

def blast_best_hit_list(blast, tools, cutoff):
    if not os.path.exists(blast):
        print '{} not exists'.format(blast)
        sys.exit(2)
    m8 = 0
    if blast.endswith('.m8'):
        m8 = 1

    solar_out = '{}.solar'.format(os.path.splitext(blast)[0])
    best_id   = '{}.solar.id'.format(os.path.splitext(blast)[0])
    cmds = []
    if not os.path.exists(solar_out):
        if m8 == 1:
            cmds.append('{} -f m8 -d -1 {} > {}'.format(tools["solar"], blast, solar_out))
        else:
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
    gene_fusion_list = 'optimize_gene_model.gene_fusion.manual.list'
 
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
    gene_low_quality_id, gene_in_repeat_id = filter_repeat_gene(args.evm_gff, args.repeat_gff, gene_swissprot_id, gene_exon_num)
    evm_noTE_highqual_gff  = write_gene_gff(evm_noTE_gff, gene_low_quality_id, args.output, 'noTE_highqual', 1)
    write_pep_from_gff(evm_noTE_highqual_gff, args.genome, tools)

    #process AS
    gene_evm_noTE_highqual_id = get_gene_id(evm_noTE_highqual_gff)
    #evm_noTE_highqual_AS_gff  = write_gene_gff(args.pasa_as_gff, gene_evm_noTE_highqual_id, args.output, 'noTE_highqual_AS', 0)
    #write_pep_from_gff(evm_noTE_highqual_AS_gff, args.genome, tools)
    evm_noTE_highqual_AS_best_gff = select_gene_model(gene_evm_noTE_highqual_id, evm_noTE_highqual_gff, args.pasa_as_gff, args.maker_gff, args.output, 'noTE_highqual_AS_best', gene_fusion_list, gene_in_repeat_id)
    write_pep_from_gff(evm_noTE_highqual_AS_best_gff, args.genome, tools) 

if __name__ == '__main__':
    main()

