#!/bin/bash
date

echo "Starting the MuonGun job"
echo "Argument line : " $@

echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

i3env=/home/hignight/work/oscNext_official/oscNext/build_trunk_jan21_py2_v3.1.1/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/src/I3FileCode/simAnalysis/fluxData/getFluxData.py
echo "Will use script: " $script

INFILE=$1
echo "Input File: " $INFILE

$i3env python $script -i $INFILE
date 
