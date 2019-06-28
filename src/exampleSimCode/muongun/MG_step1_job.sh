#!/bin/bash
date

echo "Starting the MuonGun job"
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
script=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/muongun/step_1_muongun.py
echo "Will use script: " $script

RUNTYPE=$1
ENERGYRANGE=$2
NUMEVENTS=$3
echo "Run type: " $RUNTYPE
echo "Energy range: " $ENERGYRANGE
echo "Number of events: " $NUMEVENTS

POWERLAWINDEX=-3.0
POWERLAWOFFSET=150.0
KDE=false
if [ "$RUNTYPE" == "testString" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=$4
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/testString/${GCDNAME}.i3.gz
    INCYLINDER=true 
    INCYLINDERX=0
    INCYLINDERY=0
    INCYLINDERZ=800
    INCYLINDERLENGTH=1000
    INCYLINDERRADIUS=100
    OUTCYLINDERX=0
    OUTCYLINDERY=0
    OUTCYLINDERZ=800
    OUTCYLINDERLENGTH=2000
    OUTCYLINDERRADIUS=500
elif [ "$RUNTYPE" == "HorizGeo" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=$4
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/testString/${GCDNAME}.i3.gz
    INCYLINDER=true 
    INCYLINDERX=0
    INCYLINDERY=0
    INCYLINDERZ=-600
    INCYLINDERLENGTH=200
    INCYLINDERRADIUS=800
    OUTCYLINDERX=0
    OUTCYLINDERY=0
    OUTCYLINDERZ=-600
    OUTCYLINDERLENGTH=400
    OUTCYLINDERRADIUS=1700
elif [ "$RUNTYPE" == "IceCube" ]; then
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
    echo "Using IceCube GCD"
    INCYLINDER=true 
    InCYLINDERX=46.29
    INCYLINDERY=-34.88
    INCYLINDERZ=-300.0
    INCYLINDERLENGHT=500.0
    INCYLINDERRADIUS=150.0
    OUTCYLINDERX=0
    OUTCYLINDERY=0
    OUTCYLINDERZ=0
    OUTCYLINDERLENGTH=1600
    OUTCYLINDERRADIUS=800
else 
    echo "No configuration for " $RUNTYPE "... exiting"
    exit
fi

if [ "$ENERGYRANGE" == "A" ]; then
    echo "Found configuration for " $ENERGYRANGE
    MINENERGY=1
    MAXENERGY=100
elif [ "$ENERGYRANGE" == "B" ]; then
    echo "Found configuration for " $ENERGYRANGE
    MINENERGY=100
    MAXENERGY=500
elif [ "$ENERGYRANGE" == "C" ]; then
    echo "Found configuration for " $ENERGYRANGE
    MINENERGY=500
    MAXENERGY=1000
elif [ "$ENERGYRANGE" == "D" ]; then
    echo "Found configuration for " $ENERGYRANGE
    MINENERGY=1000
    MAXENERGY=10000
elif [ "$ENERGYRANGE" == "E" ]; then
    echo "Found configuration for " $ENERGYRANGE
    MINENERGY=10000
    MAXENERGY=100000
else 
    echo "No configuration for " $ENERGYRANGE "... exiting"
    exit
fi


if [ "$INCYLINDER" = "true" ]; then 
    echo "Found settings for INCYLINDER" 
    INNERCYLINDERSETTINGS="--inner-cylinder --inner-cylinder-x "$INCYLINDERX" --inner-cylinder-y "$INCYLINDERY" --inner-cylinder-z "$INCYLINDERZ" --inner-cylinder-radius "$INCYLINDERRADIUS" --inner-cylinder-length "$INCYLINDERLENGTH
else
    echo "Will not use cylinder" 
    INNERCYLINDERSETTINGS=""
fi

OUTERCYLINDERSETTINGS="--outer-cylinder-x "$OUTCYLINDERX" --outer-cylinder-y "$OUTCYLINDERY" --outer-cylinder-z "$OUTCYLINDERZ" --outer-cylinder-radius "$OUTCYLINDERRADIUS" --outer-cylinder-length "$OUTCYLINDERLENGTH

if [ "$KDE" = "true" ]; then 
    echo "Turning on KDE " 
    KDESETTING=" --kde"
else
    echo "No KDE will be used" 
    KDESETTING=""
fi

echo "MINENERGY        : "$MINENERGY
echo "MAXENERGY        : "$MAXENERGY
echo "POWER LAW INDEX  : "$POWERLAWINDEX
echo "NUMBER OF EVENTS : "$NUMEVENTS
echo "POWERLAW OFFSET  : "$POWERLAWOFFSET
echo "INCYLINDER LINE    : ""\""$INNERCYLINDERSETTINGS"\""
echo "OUTCYLINDER LINE    : ""\""$OUTERCYLINDERSETTINGS"\""\

OUTNAME=MuonGun_step1_${RUNTYPE}_${FILE_NR}.i3.bz2
OUTDIR=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/I3Files/muongun/muongun_step1
echo "OUTFILE NAME : " ${OUTNAME}
$i3env python $script -o ${OUTDIR}/${OUTNAME} -g ${GCD_FILE} -d ${RUNTYPE} --f ${FILE_NR} --numevents ${NUMEVENTS} --min-energy ${MINENERGY} --max-energy ${MAXENERGY} --power-law-index ${POWERLAWINDEX} --power-law-offset ${POWERLAWOFFSET} ${INNERCYLINDERSETTINGS} ${KDESETTING} ${OUTERCYLINDERSETTINGS}

date 