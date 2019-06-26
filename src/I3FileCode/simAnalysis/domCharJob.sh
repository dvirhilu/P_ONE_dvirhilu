#!/bin/bash
date

echo "Starting the domChar job"
echo "Argument line : " $@

# source /home/hignight/setup_oscnext.sh
echo "TASK ID " $SLURM_ARRAY_TASK_ID
FILE_NR=`expr $SLURM_ARRAY_TASK_ID - 1`
FILE_NR=`printf "%06d\n" $FILE_NR`
echo "Filename ID : " $FILE_NR

echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

i3env=/home/hignight/work/oscNext_official/oscNext/build_trunk_jan21_py2_v3.1.1/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/src/I3FileCode/simAnalysis/createDOMCharacteristicsFiles.py
echo "Will use script: " $script

echo "OUTFILE NAME : " ${OUTNAME}
$i3env python $script
date 