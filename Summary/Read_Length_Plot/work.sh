echo "testing"
perl CompareReadLength.py --fasta ecoli/preads4falcon.fasta --list ecoli/input.fofn --project Pacbio_FALCON

echo "Citrus first 8 SMARTcell"
perl CompareReadLength.py --fasta Citrus/Citrus.fasta --list Citrus/input.fofn --project Pacbio_Citrus_8

echo "Cistus second 16 SMARTcell"
perl CompareReadLength.py --list Citrus/input.16.fofn --project Pacbio_Citrus_16
python SummaryLibReadLength.py --list Citrus/input.16.fofn --project Pacbio_Citrus_16 > log 2>&1 &
python SummaryLibReadLength.py --list Citrus/input.fofn --project Pacbio_Citrus_8 > log 2>&1 &

