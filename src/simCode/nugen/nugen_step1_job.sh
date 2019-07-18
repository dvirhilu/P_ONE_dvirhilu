#!/bin/bash
date
startsecond=$(date +%s)

echo "Starting the NuGen job"
echo "Argument line : " $@

# source /home/hignight/setup_oscnext.sh
#echo "TASK ID " $SLURM_ARRAY_TASK_ID
#FILE_NR=`expr $SLURM_ARRAY_TASK_ID - 1`
#FILE_NR=`printf "%06d\n" $FILE_NR`
FILE_NR=$6
echo "Filename ID : " $FILE_NR

echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

i3env=/home/users/hignight/oscnext/build_trunk_july_02_2019/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/users/dhilu/P_ONE_dvirhilu/src/simCode/nugen/step1_neutrino_generator.py
echo "Will use script: " $script

RUNTYPE=$1
NUMEVENTS=$2
LOGMINENERGY=$3
LOGMAXENERGY=$4
CONENAGLE=$5

echo "Run type: " $RUNTYPE
echo "Energy range: " $ENERGYRANGE
echo "Number of events: " $NUMEVENTS

if [ "$RUNTYPE" == "testString" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=HorizTestString_n15_b100.0_v50.0_l1_simple_spacing.i3.gz
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/home/users/dhilu/I3Files/gcd/testStrings/${GCDNAME}.i3.gz
 
    CYLINDERX=0
    CYLINDERY=0
    CYLINDERZ=116.08
    CYLINDERLENGTH=1600
    CYLINDERRADIUS=800
elif [ "$RUNTYPE" == "HorizGeo" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=HorizGeo_n10_b100.0_a90.0_l1_linear_reset_offset_exp_r_spacing
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/home/users/dhilu/I3Files/gcd/uncorHorizGeo/${GCDNAME}.i3.gz

    CYLINDERX=0
    CYLINDERY=0
    CYLINDERZ=-600
    CYLINDERLENGTH=400
    CYLINDERRADIUS=1700
elif [ "$RUNTYPE" == "IceCube" ]; then
    echo "Found configuration for " $RUNTYPE
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/home/users/dhilu/I3Files/gcd/IceCube/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
    echo "Using IceCube GCD"

    CYLINDERX=0
    CYLINDERY=0
    CYLINDERZ=0
    CYLINDERLENGTH=1600
    CYLINDERRADIUS=800
elif [ "$RUNTYPE" == "cube" ]; then
    echo "Found configuration for " $RUNTYPE
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/home/users/dhilu/I3Files/gcd/cube/cubeGeometry_1600_15_50.i3.gz
    echo "Using cube GCD"
    
    CYLINDERX=0
    CYLINDERY=0
    CYLINDERZ=0
    CYLINDERLENGTH=1600
    CYLINDERRADIUS=800
else 
    echo "No configuration for " $RUNTYPE "... exiting"
    exit
fi

CYLINDERSETTINGS="-x "$CYLINDERX" -y "$CYLINDERY" -z "$CYLINDERZ" -r "$CYLINDERRADIUS" -l "$CYLINDERLENGTH
POWERLAWINDEX=2
OUTNAME=NuGen_step1_${RUNTYPE}_${FILE_NR}.i3.gz
OUTDIR=/home/users/dhilu/I3Files/nugenStep1

echo "FILE NUMBER      : "$FILE_NR
echo "NUMBER OF EVENTS : "$NUMEVENTS
echo "OUTPUT FILE NAME : "$OUTNAME
echo "OUTPUT FILE DIR  : "$OUTDIR
echo "Keeping flavours, ratios, zenithRange as default"
echo "LOG ENERGY RANGE : "${LOGMINENERGY}":"${LOGMAXENERGY}
echo "POWER LAW INDEX  : "$POWERLAWINDEX
echo "GENERATING EVENTS IN A CONE ANGLE: "$CONEANGLE

echo "CYLINDER LINE    : ""\""$CYLINDERSETTINGS"\""

$i3env python $script -N ${FILE_NR} -n $NUMEVENTS -o ${OUTDIR}/${OUTNAME} -s ${FILE_NR}000 -g $GCD_FILE -a $CONEANGLE -E ${LOGMINENERGY}":"${LOGMAXENERGY} -p ${POWERLAWINDEX} ${CYLINDERSETTINGS}

date 
endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
