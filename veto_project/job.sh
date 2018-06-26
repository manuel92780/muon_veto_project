#!/bin/bash

#prep inputs
num_events=$1
event_start=$2
event_end=$3
num=$(($event_start / $num_events))

#paths
dir=/home/msilva/muon_gun_multiplicity_study/veto_project
mkdir $dir/data
scratch_dir=/home/msilva/scratch
mkdir $scratch_dir
mkdir $scratch_dir/muon_gun_veto
cd $scratch_dir/muon_gun_veto/
cp $dir/generate_muons.py .

#setup env
eval `/cvmfs/icecube.opensciencegrid.org/py2-v2/setup.sh`

#run job
bash /home/msilva/combo/build/env-shell.sh <<EOF
python generate_muons.py --nseed ${num} --use-gpu --nevents $num_events --out $dir/data/muons_${event_start}_${event_end}.i3.gz
EOF

cd $dir

