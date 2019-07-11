#!/bin/bash
date

echo "Starting the MuonGun job"
echo "Argument line : " $@

#echo "TASK ID " $SLURM_ARRAY_TASK_ID
#FILE_NR=`expr $SLURM_ARRAY_TASK_ID - 1`
#FILE_NR=`printf "%06d\n" $FILE_NR`
#echo "Filename ID : " $FILE_NR
FILE_NR=000900
ARR_NUM=901

echo ${OPENCL_VENDOR_PATH}
echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo ${OPENCL_VENDOR_PATH}

i3env=/home/hignight/work/oscNext_official/oscNext/build_trunk_jan21_py2_v3.1.1/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/project/6008051/dvirhilu/P_ONE_dvirhilu/src/simCode/nugen/step2_clsim_setCrossE.py
echo "Will use script: " $script

RUNTYPE=$1
MEDIUMNAME=$2

if [ "$RUNTYPE" == "testString" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=TestString_n10_b50.0_v50.0_l1_simple_spacing
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/testStrings/${GCDNAME}.i3.gz
 
elif [ "$RUNTYPE" == "HorizGeo" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=HorizGeo_n10_b100.0_a90.0_l1_linear_reset_offset_exp_r_spacing
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/uncorHorizGeo/${GCDNAME}.i3.gz
 
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

MEDIUMMODEL=/project/6008051/dvirhilu/P_ONE_dvirhilu/propagationMediumModels/${MEDIUMNAME}
CROSSENERGY=200.0
echo "MEDIUMMODEL    : " $MEDIUMMODEL
echo "CROSSENERGY : " $CROSSENERGY
#INFILENAME=NuGen_step1_${RUNTYPE}_${FILE_NR}.i3.gz
INFILENAME=testFile.i3.gz
INFOLDER=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/I3Files/nugen/nugenStep1

echo "INFILEPATH: " ${INFOLDER}/$INFILENAME
OUTFILENAME=NuGen_step2_${RUNTYPE}_test.i3.gz
OUTFOLDER=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/I3Files/nugen/nugenStep2
echo "OUTFILEPATH : " ${OUTFOLDER}/$OUTFILENAME

echo "GCD: " $GCD_FILE

$i3env python $script -t -i ${INFOLDER}/$INFILENAME -g $GCD_FILE -o ${OUTFOLDER}/$OUTFILENAME -r 14000 -m $MEDIUMMODEL -c $CROSSENERGY -l $ARR_NUM
#$i3env python $script -i ${INFOLDER}/$INFILENAME -g $GCD_FILE -o ${OUTFOLDER}/$OUTFILENAME -r 14000 -m $MEDIUMMODEL -c $CROSSENERGY -l $ARR_NUM
