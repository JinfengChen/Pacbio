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
python EVM_from_Maker.py --gff Fairchildv1.all.gff --genome Fairchildv1.fa

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

def maker(gffile, fastafile, repeat, outdir):
    """
    %prog maker maker.gff3 genome.fasta

    Prepare EVM inputs by separating tracks from MAKER.
    """
    from jcvi.formats.base import SetFile, FileShredder
    from jcvi.formats.base import write_file
    from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh

    A, T, P = "ABINITIO_PREDICTION", "TRANSCRIPT", "PROTEIN"
    # Stores default weights and types
    Registry = {\
        "maker": (A, 5),
        "augustus_masked": (A, 1),
        "snap_masked": (A, 1),
        "genemark": (A, 1),
        "fgenesh_masked": (A, 1),
        "est2genome": (T, 5),
        "est_gff": (T, 5),
        "protein2genome": (P, 5),
        "blastx": (P, 1)
    }

    types = "type.ids"
    if need_update(gffile, types):
        cmd = "cut -f2 -s {0} | sort -u".format(gffile)
        sh(cmd, outfile=types)

    types = SetFile(types)
    reg = defaultdict(list)
    weightsfile = "weight.txt"
    contents = []
    for s in types:
        rs = s.split(":")[0]
        if rs not in Registry:
            continue

        type, weight = Registry[rs]
        reg[type].append(s)
        contents.append("\t".join(str(x) for x in (type, s, weight)))

    contents = "\n".join(sorted(contents))
    write_file(weightsfile, contents)

    evs = [x + ".gff" for x in (A, T, P)]
    FileShredder(evs)

    for type, tracks in reg.items():
        for t in tracks:
            if type == T or type == P:
                cmd = "grep '\t{0}' {1} | grep -v '_match\t' | sed 's/ID=.*;Parent=/ID=/' >> {2}.gff".format(t, gffile, type)
                sh(cmd)
            elif t == "maker":
                cmd = "grep '\t{0}' {1} | grep -v '_match\t' >> {2}.gff".format(t, gffile, type)
                sh(cmd)
            else:
                cmd = "grep '\t{0}' {1} | grep -v '_match\t' >> {2}.temp.gff".format(t, gffile, type)
                sh(cmd)
 
    predict_gff(A)
    sh('mkdir {}'.format(outdir))
    sh('mv {}.* {}'.format(A, outdir))
    sh('mv {}.* {}'.format(P, outdir))
    sh('mv {}.* {}'.format(T, outdir))
    sh('mv {}.* {}'.format('weight.txt', outdir))
    sh('ln -s {} {}/'.format(os.path.abspath(fastafile), outdir)) 
    sh('ln -s {} {}/'.format(os.path.abspath(repeat), outdir))   


'''
scaffold74      genemark        match   8389    10530   .       +       .       ID=scaffold74:hit:45659:4.5.0.0;Name=maker-scaffold74-pred_gff_genemark-gene-0.8-mRNA-1
scaffold74      genemark        match_part      8389    8541    .       +       .       ID=scaffold74:hsp:69250:4.5.0.0;Parent=scaffold74:hit:45659:4.5.0.0;Target=maker-scaffold74-pred_gff_genemark-gene-0.8-mRNA-1 1 153 +;Gap=M153
scaffold74      genemark        match_part      8681    8827    .       +       .       ID=scaffold74:hsp:69251:4.5.0.0;Parent=scaffold74:hit:45659:4.5.0.0;Target=maker-scaffold74-pred_gff_genemark-gene-0.8-mRNA-1 154 300 +;Gap=M
scaffold74      genemark        match_part      9417    10064   .       +       .       ID=scaffold74:hsp:69252:4.5.0.0;Parent=scaffold74:hit:45659:4.5.0.0;Target=maker-scaffold74-pred_gff_genemark-gene-0.8-mRNA-1 301 948 +;Gap=M
scaffold74      genemark        match_part      10381   10530   .       +       .       ID=scaffold74:hsp:69253:4.5.0.0;Parent=scaffold74:hit:45659:4.5.0.0;Target=maker-scaffold74
'''
def predict_gff(A):
    gff_temp = '{}.temp.gff'.format(A)
    gff      = '{}.gff'.format(A)
    ofile = open(gff, 'a')
    with open (gff_temp, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'): 
                unit  = re.split(r'\t',line)
                tool    = unit[1]
                feature = unit[2]
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    #print attr
                    if not attr == '':
                        #print 'yes'
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                if feature == "match":
                    gid   = temp["ID"]
                    gname = temp["Name"]
                    #gene
                    unit[2] = 'gene'
                    unit[8] = 'ID={}.gene;Name={}'.format(gid, gname)
                    print >> ofile, '\t'.join(unit)
                    #mRNA
                    unit[2] = 'mRNA'
                    unit[8] = 'ID={}.mRNA;Parent={}.gene'.format(gid, gid)
                    print >> ofile, '\t'.join(unit)
                elif feature == "match_part":
                    eid   = temp["ID"]
                    gid   = temp["Parent"] 
                    #exon
                    unit[2] = 'exon'
                    unit[8] = 'ID={}:exon;Parent={}.mRNA'.format(eid, gid)
                    print >> ofile, '\t'.join(unit)
                    #CDS
                    unit[2] = 'CDS'
                    unit[8] = 'ID={}:cds;Parent={}.mRNA'.format(eid, gid)
                    print >> ofile, '\t'.join(unit)
            else:
                print >> ofile, line
    ofile.close()

def write_slurm(genome, repeat, outdir):
    cmd='''
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=100G
#SBATCH --time=40:00:00
#SBATCH --output=EVM.sh.stdout
#SBATCH -p intel
#SBATCH --workdir=./


start=`date +%%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

genome=%s
repeat=%s
predict=ABINITIO_PREDICTION.gff
est2genome=TRANSCRIPT.gff
protein2genome=PROTEIN.gff

module load EVM/1.1.1
#EVM_dir=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Annotation/gene/EVM/1.1.1
EVM_dir=/bigdata/bioinfo/pkgadmin/opt/linux/centos/7.x/x86_64/pkgs/EVM/1.1.1/
#perl $EVM_dir/evidence_modeler.pl --weights weight.txt --genome $genome --gene_predictions $predict --protein_alignments $protein2genome --transcript_alignments $est2genome --repeat $repeat
#perl /opt/linux/centos/7.x/x86_64/pkgs/EVM/1.1.1/EvmUtils/EVM_to_GFF3.pl evm.out.orig chr1 > evm.out.gff3

#partition
perl $EVM_dir/EvmUtils/partition_EVM_inputs.pl --genome $genome --gene_predictions $predict --protein_alignments $protein2genome --transcript_alignments $est2genome --repeat $repeat --segmentSize 200000 --overlapSize 20000 --partition_listing partitions_list.out

perl $EVM_dir/EvmUtils/write_EVM_commands.pl --genome $genome --gene_predictions $predict --protein_alignments $protein2genome --transcript_alignments $est2genome --repeat $repeat --weights `pwd`/weight.txt --output_file_name evm.out  --partitions partitions_list.out > commands.list

perl /rhome/cjinfeng/BigData/software/bin/qsub-slurm.pl --maxjob 90 --lines 10 --interval 120 --task 1 --mem 30G --time 20:00:00 --convert no commands.list

perl $EVM_dir/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

perl $EVM_dir/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome $genome

find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3


end=`date +%%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

''' %(genome, repeat)

    ofile = open('{}/run_EVM.sh'.format(outdir), 'w')
    print >> ofile, cmd
    ofile.close()
     
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff')
    parser.add_argument('--genome')
    parser.add_argument('--repeat')
    parser.add_argument('--outdir')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        os.path.exists(args.gff) and os.path.exists(args.genome)
    except:
        usage()
        sys.exit(2)

    if not args.outdir:
        args.outdir = 'run_EVM_FCM'

    #maker(args.gff, args.genome, args.repeat, args.outdir)
    write_slurm(args.genome, args.repeat, args.outdir)

if __name__ == '__main__':
    main()

