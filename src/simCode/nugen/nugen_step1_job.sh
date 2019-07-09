#!/bin/bash
date

echo "Starting the MuonGun job"
echo "Argument line : " $@

# source /home/hignight/setup_oscnext.sh
#echo "TASK ID " $SLURM_ARRAY_TASK_ID
#FILE_NR=`expr $SLURM_ARRAY_TASK_ID - 1`
#FILE_NR=`printf "%06d\n" $FILE_NR`
FILE_NR=900
echo "Filename ID : " $FILE_NR

echo "Starting cvmfs " 
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

i3env=/home/hignight/work/oscNext_official/oscNext/build_trunk_jan21_py2_v3.1.1/env-shell.sh
echo "Will use i3 environment: " ${i3env}
script=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/src/simCode/nugen/step1_neutrino_generator.py
echo "Will use script: " $script

RUNTYPE=$1
NUMEVENTS=$2
LOGMINENERGY=$3
LOGMAXENERGY=$4

echo "Run type: " $RUNTYPE
echo "Energy range: " $ENERGYRANGE
echo "Number of events: " $NUMEVENTS

if [ "$RUNTYPE" == "testString" ]; then
    echo "Found configuration for " $RUNTYPE
    CYLINDERX=0
    CYLINDERY=0
    CYLINDERZ=0
    CYLINDERLENGTH=1600
    CYLINDERRADIUS=300
elif [ "$RUNTYPE" == "HorizGeo" ]; then
    echo "Found configuration for " $RUNTYPE
    CYLINDERX=0
    CYLINDERY=0
    CYLINDERZ=-600
    CYLINDERLENGTH=400
    CYLINDERRADIUS=1700
elif [ "$RUNTYPE" == "IceCube" ]; then
    echo "Found configuration for " $RUNTYPE
    CYLINDERX=0
    CYLINDERY=0
    CYLINDERZ=0
    CYLINDERLENGTH=1600
    CYLINDERRADIUS=800
elif [ "$RUNTYPE" == "cube" ]; then
    echo "Found configuration for " $RUNTYPE
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
POWERLAWINDEX=2.0
OUTNAME=NuGen_step1_${RUNTYPE}_${FILE_NR}.i3.gz
OUTDIR=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/I3Files/nugen/nugenStep1/

echo "FILE NUMBER      : "$FILE_NR
echo "NUMBER OF EVENTS : "$NUMEVENTS
echo "OUTPUT FILE NAME : "$OUTNAME
echo "OUTPUT FILE DIR  : "$OUTDIR
echo "Keeping flavours, ratios, zenithRange as default"
echo "LOG ENERGY RANGE : "${LOGMINENERGY}":"${LOGMAXENERGY}
echo "POWER LAW INDEX  : "$POWERLAWINDEX

echo "CYLINDER LINE    : ""\""$CYLINDERSETTINGS"\""

$i3env python $script -N ${FILE_NR} -n $NUMEVENTS -o ${OUTDIR}/${OUTNAME} -s ${FILE_NR}000 -E ${LOGMINENERGY}":"${LOGMAXENERGY} -p ${POWERLAWINDEX} ${CYLINDERSETTINGS}

date 
