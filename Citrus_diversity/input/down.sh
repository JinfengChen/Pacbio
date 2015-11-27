#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

#Seville sour orange WGS SSO SRX372786
#wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR102/SRR1023728/SRR1023728.sra
#mv SRR1023728.sra SSO.sra
#Low acid pummelo WGS LAP SRX372702
#/opt/aspera/3.3.3/bin/ascp -i /opt/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -k 1 -T -l40m anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR102/SRR1023638/SRR1023638.sra ./
#mv SRR1023638.sra LAP.sra
#W. Murcott mandarin WGS WMM SRX372687
#/opt/aspera/3.3.3/bin/ascp -i /opt/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -k 1 -T -l40m anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR102/SRR1023626/SRR1023626.sra ./
#mv SRR1023626.sra WMM.sra
#Sweet orange WGS SWO SRX372703
#/opt/aspera/3.3.3/bin/ascp -i /opt/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -k 1 -T -l40m anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR102/SRR1023639/SRR1023639.sra ./
#mv SRR1023639.sra SWO.sra
#Clementine mandarin CLM WGS SRX371962
#/opt/aspera/3.3.3/bin/ascp -i /opt/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -k 1 -T -l40m anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR102/SRR1022654/SRR1022654.sra ./
#mv SRR1022654.sra CLM.sra
#Willowleaf mandarin WGS WLM SRX372685
/opt/aspera/3.3.3/bin/ascp -i /opt/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -k 1 -T -l40m anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR102/SRR1023625/SRR1023625.sra ./
mv SRR1023625.sra WLM.sra
#Chandler pummelo WGS CHP SRX372688
/opt/aspera/3.3.3/bin/ascp -i /opt/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -k 1 -T -l40m anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR102/SRR1023627/SRR1023627.sra ./
mv SRR1023627.sra CHP.sra
#Ponkan mandarin WGS PKM SRX372665
#/opt/aspera/3.3.3/bin/ascp -i /opt/aspera/3.3.3/etc/asperaweb_id_dsa.openssh -k 1 -T -l20m anonftp@ftp.ncbi.nlm.nih.gov:sra/sra-instant/reads/ByRun/sra/SRR/SRR102/SRR1023619/SRR1023619.sra ./
#mv SRR1023619.sra PKM.sra

echo "Done"
