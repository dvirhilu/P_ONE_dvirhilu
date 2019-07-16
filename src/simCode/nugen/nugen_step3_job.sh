#!/bin/bash
date

echo "Starting the NuGen job"
echo "Argument line : " $@

startsecond=$(date +%s)

echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

MIN_FILE_NR=$1
MAX_FILE_NR=$2
HIT_NUM_THRESH=$3
DOM_NUM_THRESH=$4
GCDTYPE=$5
SIMTYPE=nugen
DOMTYPE=MDOM

echo "Runs will go from file number: " $MIN_FILE_NR " to " $MAX_FILE_NR
echo "Number of hits needed to preserve frame: "$HIT_NUM_THRESH
echo "Number of DOMs with hits needed to preserve frame: " $DOM_NUM_THRESH
echo "Using GCD for: "$GCDTYPE
echo "Simulations done with: "$SIMTYPE

for value in $(seq $MIN_FILE_NR $MAX_FILE_NR); do
    python generateHitsFromI3Photons.py -n $value -s $SIMTYPE -g $GCDTYPE -d $DOMTYPE -H $HIT_NUM_THRESH -D $DOM_NUM_THRESH
    echo "finished loop for: "$value
done


date 
endsecond=$(date +%s)
echo "This job took : "`expr $endsecond - $startsecond`" s"
