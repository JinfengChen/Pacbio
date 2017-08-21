from pybedtools import BedTool
repeats = BedTool('Fairchild_v1.fasta.RepeatMasker.out.gff') 
genes   = BedTool('Fairchild.optimized_model.noTE_highqual_AS.gff')
for g in genes.intersect(repeats):
    print g
