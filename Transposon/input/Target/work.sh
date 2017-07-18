
echo "Run Target on MITEhunter results"
python make_target_nonauto_general_redo.py /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Transposon/input/Target/query/ /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Transposon/input/Target/Target /rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/Transposon/input/Target/reference/Fairchild_hapltype1.fasta Target_Run
echo "Quick summary of high identical families"
python Select_Multiple_HighIdentity.py --input Target/Target_Run_2017_02_15_132003
