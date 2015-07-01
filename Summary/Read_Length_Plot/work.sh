echo "testing"
perl CompareReadLength.py --fasta ecoli/preads4falcon.fasta --list ecoli/input.fofn --project Pacbio_FALCON

echo "Citrus first 8 SMARTcell"
perl CompareReadLength.py --fasta Citrus/Citrus.fasta --list Citrus/input.fofn --project Pacbio_Citrus_8

echo "Cistus second 16 SMARTcell"
perl CompareReadLength.py --list Citrus/input.16.fofn --project Pacbio_Citrus_16
python SummaryLibReadLength.py --list Citrus/input.16.fofn --project Pacbio_Citrus_16 > log 2>&1 &
python SummaryLibReadLength.py --list Citrus/input.fofn --project Pacbio_Citrus_8 > log 2>&1 &

echo "Citrus third 12 SMRT cell"
perl CompareReadLength.py --list Citrus/input.12.fofn --project Pacbio_Citrus_12
python SummaryLibReadLength.py --list Citrus/input.12.fofn --project Pacbio_Citrus_12 > log 2>&1 &

echo "Citrus forth 4 SMRT cell"
perl CompareReadLength.py --list Citrus/input.4.fofn --project Pacbio_Citrus_4
python SummaryLibReadLength.py --list Citrus/input.4.fofn --project Pacbio_Citrus_4 > log 2>&1 &

echo "Corrected reads 36 SMRT cell"
perl CompareReadLength.py --fasta Citrus/preads4falcon.fasta --list Citrus/input.36.fofn --project Pacbio_Citrus_36
perl CompareReadLength.py --list Citrus/input_preads.fofn --project Pacbio_Citrus_36_preads

echo "Corrected read 40 SMRT cell"
perl CompareReadLength.py --fasta Citrus/preads4falcon.fasta --list Citrus/input.40.fofn --project Pacbio_Citrus_40

