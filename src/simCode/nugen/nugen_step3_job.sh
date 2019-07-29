#!/bin/bash
date

echo "Starting the NuGen job"
echo "Argument line : " $@

startsecond=$(date +%s)

echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`



RUNNUM=$1
HIT_NUM_THRESH=$2
DOM_NUM_THRESH=$3
GCDTYPE=$4
ISLOCAL=$5
SIMTYPE=nugen
DOMTYPE=MDOM

if [ "$ISLOCAL" == "true" ]; then
    echo "Running locally"
    FILEPATH=/home/dvir/workfolder/
    echo "File path: "$FILEPATH

    i3env=/home/dvir/combo/build/env-shell.sh
    echo "Will use i3 environment: " ${i3env}
 
elif [ "$ISLOCAL" == "false" ]; then
    echo "Running on illume"
    FILEPATH=/home/users/dhilu/
    echo "File path: "$FILEPATH

    i3env=/home/users/hignight/oscnext/build_trunk_july_02_2019/env-shell.sh
    echo "Will use i3 environment: " ${i3env}
 
else 
    echo "No configuration for " $ISLOCAL "... exiting"
    exit
fi

script=${FILEPATH}P_ONE_dvirhilu/src/simCode/nugen/generateHitsFromI3Photons.py
echo "Will use script: " $script

echo "Run number: " $RUNNUM
echo "Number of hits needed to preserve DOM: "$HIT_NUM_THRESH
echo "Number of DOMs with hits needed to preserve frame: " $DOM_NUM_THRESH
echo "Using GCD for: "$GCDTYPE
echo "Simulations done with: "$SIMTYPE

$i3env python $script -n $RUNNUM -s $SIMTYPE -g $GCDTYPE -d $DOMTYPE -H $HIT_NUM_THRESH -D $DOM_NUM_THRESH -f $FILEPATH

date 
endsecond=$(date +%s)
echo "This job took : "`expr $endsecond - $startsecond`" s"
