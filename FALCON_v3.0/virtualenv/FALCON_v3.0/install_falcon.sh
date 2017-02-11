echo "FALCON"

unset PYTHONPATH # must run before runing fc_run
unset FALCON_WORKSPACE
unset PYTHONUSERBASE
unset FALCON_PREFIX

git clone git://github.com/PacificBiosciences/FALCON-integrate.git
cd FALCON-integrate
git checkout master  # or whatever version you want
make init
source env.sh
make config-edit-user
bash ./FALCON-make/config-edit-user.sh
make -j all
make test  # to run a simple one

echo "FALCON update"
unset PYTHONPATH
unset FALCON_WORKSPACE
unset PYTHONUSERBASE
unset FALCON_PREFIX
make init
git submodule update --init
cp -f default-env.sh env.sh
source env.sh
make config-edit-user
bash ./FALCON-make/config-edit-user.sh
make -j all
make test

echo "FALCON_unzip"
git clone https://github.com/PacificBiosciences/FALCON_unzip.git
cd FALCON_unzip
export PYTHONPATH=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON_v3.0/virtualenv/FALCON_v3.0/FALCON-integrate/fc_env//lib/python2.7/site-packages/
python setup.py install --prefix ~/BigData/00.RD/Assembly/Pacbio/FALCON_v3.0/virtualenv/FALCON_v3.0/FALCON-integrate/fc_env/

echo "pysam"
export PYTHONPATH=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON_v3.0/virtualenv/FALCON_v3.0/FALCON-integrate/fc_env/lib64/python2.7/site-packages
python setup.py install --prefix=/rhome/cjinfeng/BigData/00.RD/Assembly/Pacbio/FALCON_v3.0/virtualenv/FALCON_v3.0/FALCON-integrate/fc_env/
