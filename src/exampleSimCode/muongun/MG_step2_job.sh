#!/bin/bash
date

echo "Starting the MuonGun job"
echo "Argument line : " $@

echo "TASK ID " $SLURM_ARRAY_TASK_ID
FILE_NR=`expr $SLURM_ARRAY_TASK_ID - 1`
FILE_NR=`printf "%06d\n" $FILE_NR`
echo "Filename ID : " $FILE_NR

echo ${OPENCL_VENDOR_PATH}
echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo ${OPENCL_VENDOR_PATH}

i3env=/home/hignight/work/oscNext_official/oscNext/build_trunk_jan21_py2_v3.1.1/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/projects/6008051/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/muongun/step_2_clsim_setCrossE.py
echo "Will use script: " $script

OUTRUN=$1

if [ "$OUTRUN" == "139005" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139005
    ICEMODEL="ANTARES"
    CROSSENERGY=200.0
elif [ "$OUTRUN" == "139006" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139006
    ICEMODEL="ANTARES"
    CROSSENERGY=200.0
elif [ "$OUTRUN" == "139008" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139008
    ICEMODEL="ANTARES"
    CROSSENERGY=200.0
elif [ "$OUTRUN" == "139010" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139010
    ICEMODEL="ANTARES"
    CROSSENERGY=200.0
else 
    echo "No configuration for " $OUTRUN "... exiting"
    exit
fi
echo "INRUN       : " $INRUN
echo "ICEMODEL    : " $ICEMODEL
echo "CROSSENERGY : " $CROSSENERGY
INFILENAME=MuonGun_step1_${INRUN}_${FILE_NR}.i3.bz2
INFOLDER=/projects/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/muongun_step1

echo "INFILEPATH: " ${INFOLDER}/$INFILENAME
OUTFILENAME=MuonGun_step2_${OUTRUN}_${FILE_NR}.i3.bz2
OUTFOLDER = /projects/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/muongun_step2
echo "OUTFILEPATH : " ${OUTFOLDER}/$OUTFILENAME
GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
echo "GCD: " $GCD_FILE

$i3env python $script -t -i ${INFOLDER}/$INFILENAME -g $GCD_FILE -o ${OUTFOLDER}/$OUTFILENAME -r ${OUTRUN} -l $SLURM_ARRAY_TASK_ID -c $CROSSENERGY
