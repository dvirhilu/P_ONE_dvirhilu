#!/bin/bash
date
startsecond=$(date +%s)

echo "Starting the clsim job"
echo "Argument line : " $@

#echo "TASK ID " $SLURM_ARRAY_TASK_ID
#FILE_NR=`expr $SLURM_ARRAY_TASK_ID - 1`
#FILE_NR=`printf "%06d\n" $FILE_NR`
FILE_NR=$3
ARR_NUM=`expr $FILE_NR + 1`
#ARR_NUM=$SLURM_ARRAY_TASK_ID
#echo "Filename ID : " $FILE_NR

echo ${OPENCL_VENDOR_PATH}
echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo ${OPENCL_VENDOR_PATH}

i3env=/home/users/hignight/oscnext/build_trunk_july_02_2019/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/users/dhilu/P_ONE_dvirhilu/src/simCode/nugen/step2_clsim_setCrossE.py
echo "Will use script: " $script

RUNTYPE=$1
MEDIUMNAME=$2

if [ "$RUNTYPE" == "testString" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=HorizTestString_n15_b100.0_v50.0_l1_simple_spacing
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/home/users/dhilu/I3Files/gcd/testStrings/${GCDNAME}.i3.gz
 
elif [ "$RUNTYPE" == "HorizGeo" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=CorrHorizGeo_n15_b100.0_a18.0_l3_rise_fall_offset_simple_spacing
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/home/users/dhilu/I3Files/gcd/HorizGeo/${GCDNAME}.i3.gz
 
 elif [ "$RUNTYPE" == "denseGeo" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=denseGeo_n30_b50.0_a4.5_l7_linear_reset_offset_simple_spacing
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/home/users/dhilu/I3Files/gcd/denseGeo/${GCDNAME}.i3.gz

elif [ "$RUNTYPE" == "IceCube" ]; then
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/home/users/dhilu/I3Files/gcd/icecube/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
    echo "Using IceCube GCD"

elif [ "$RUNTYPE" == "cube" ]; then
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/home/users/dhilu/I3Files/gcd/cube/cubeGeometry_1600_15_50.i3.gz
    echo "Using cube GCD"
    
else 
    echo "No configuration for " $RUNTYPE "... exiting"
    exit
fi

MEDIUMMODEL=/home/users/dhilu/P_ONE_dvirhilu/propagationMediumModels/${MEDIUMNAME}
CROSSENERGY=200.0
echo "MEDIUMMODEL    : " $MEDIUMMODEL
echo "CROSSENERGY : " $CROSSENERGY
INFILENAME=NuGen_step1_${RUNTYPE}_${FILE_NR}.i3.gz
#INFILENAME=testFile.i3.gz
INFOLDER=/home/users/dhilu/I3Files/nugen/nugenStep1/denseGeo

echo "INFILEPATH: " ${INFOLDER}/$INFILENAME
OUTFILENAME=NuGen_step2_${RUNTYPE}_${FILE_NR}.i3.gz
OUTFOLDER=/home/users/dhilu/I3Files/nugen/nugenStep2/denseGeo
echo "OUTFILEPATH : " ${OUTFOLDER}/$OUTFILENAME

echo "GCD: " $GCD_FILE

$i3env python $script -t -i ${INFOLDER}/$INFILENAME -g $GCD_FILE -o ${OUTFOLDER}/$OUTFILENAME -r 14000 -m $MEDIUMMODEL -c $CROSSENERGY -l $ARR_NUM
#$i3env python $script -i ${INFOLDER}/$INFILENAME -g $GCD_FILE -o ${OUTFOLDER}/$OUTFILENAME -r 14000 -m $MEDIUMMODEL -c $CROSSENERGY -l $ARR_NUM

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
