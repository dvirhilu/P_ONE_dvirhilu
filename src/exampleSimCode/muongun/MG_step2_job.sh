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
script=/project/6008051/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/muongun/step_2_clsim_setCrossE.py
echo "Will use script: " $script

RUNTYPE=$1
MEDIUMNAME=$2

if [ "$RUNTYPE" == "testString" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=$3
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/testString/${GCDNAME}.i3.gz
elif [ "$RUNTYPE" == "HorizGeo" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=$3
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/testString/${GCDNAME}.i3.gz
elif [ "$RUNTYPE" == "IceCube" ]; then
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
    echo "Using IceCube GCD"
elif [ "$RUNTYPE" == "cube" ]; then
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/cube/cubeGeometry_1600_15_50.i3.gz
    echo "Using cube GCD"
else 
    echo "No configuration for " $RUNTYPE "... exiting"
    exit
fi

MEDIUMMODEL=/project/6008051/dvirhilu/P_ONE_dvirhilu/propagationMediumModels/${MEDIUMMODEL}
CROSSENERGY=200.0
echo "INRUN       : " $INRUN
echo "MEDIUMMODEL    : " $MEDIUMMODEL
echo "CROSSENERGY : " $CROSSENERGY
INFILENAME=MuonGun_step1_${RUNTYPE}_${FILE_NR}.i3.bz2
INFOLDER=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/I3Files/muongun/muongun_step1

echo "INFILEPATH: " ${INFOLDER}/$INFILENAME
OUTFILENAME=MuonGun_step2_${RUNTYPE}_${FILE_NR}.i3.bz2
OUTFOLDER=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/I3Files/muongun/muongun_step2
echo "OUTFILEPATH : " ${OUTFOLDER}/$OUTFILENAME
GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
echo "GCD: " $GCD_FILE

$i3env python $script -t -i ${INFOLDER}/${INFILENAME} -g ${GCD_FILE} -o ${OUTFOLDER}/${OUTFILENAME} -r ${RUNTYPE} -l ${SLURM_ARRAY_TASK_ID} -c ${CROSSENERGY} -m $MEDIUMMODEL
