#!/bin/bash

#prep inputs
nseed=$1
nfile=$2

#paths
dir=/data/user/msilva/corsika_20222
scratch_dir=/home/msilva/scratch
mkdir $scratch_dir
mkdir $scratch_dir/corsika
cd $scratch_dir/corsika/
cp /home/msilva/muon_gun_multiplicity_study/corsika_comparison/process_corsika.py .

#setup env
eval `/cvmfs/icecube.opensciencegrid.org/py2-v2/setup.sh`

#run job
bash /home/msilva/combo/build/env-shell.sh <<EOF
python process_corsika.py --nseed ${nseed} --out ${dir}/corsika ${nfile}  
EOF


