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

RUNNR=$1
echo "Run number: " $RUNNR

INNERCYLINDERSETTINGS=""
GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
MAXENERGY=10000.0
MINENERGY=50.0 
POWERLAWINDEX=-3.0
POWERLAWOFFSET=150.0
NUMEVENTS=1000 
KDE=false
if [ "$RUNNR" == "139005" ]; then
    echo "Found configuration for " $RUNNR
    CYLINDER=true 
    CYLINDERX=46.29
    CYLINDERY=-34.88
    CYLINDERZ=-300.0
    CYLINDERLENGHT=500.0
    CYLINDERRADIUS=150.0
#    NUMEVENTS=100000 
elif [ "$RUNNR" == "139006" ]; then
    echo "Found configuration for " $RUNNR
    CYLINDER=true 
    CYLINDERX=46.29
    CYLINDERY=-34.88
    CYLINDERZ=-300.0
    CYLINDERLENGHT=1000.0
    CYLINDERRADIUS=300.0
elif [ "$RUNNR" == "139008" ]; then
    echo "Found configuration for " $RUNNR 
    MINENERGY=50.0 
#    NUMEVENTS=860000
    CYLINDER=false 
    KDE=false
elif [ "$RUNNR" == "139010" ]; then
    echo "Found configuration for " $RUNNR 
    MINENERGY=150.0
    MAXENERGY=5000.0 
    CYLINDER=true 
    CYLINDERX=46.29
    CYLINDERY=-34.88
    CYLINDERZ=-350.0    
    CYLINDERRADIUS=200.0
    CYLINDERLENGHT=700.0
#    NUMEVENTS=250000
else 
    echo "No configuration for " $RUNNR "... exiting"
    exit
fi




if [ "$CYLINDER" = "true" ]; then 
    echo "Found settings for cylinder" 
    INNERCYLINDERSETTINGS="--inner-cylinder --inner-cylinder-x "$CYLINDERX" --inner-cylinder-y "$CYLINDERY" --inner-cylinder-z "$CYLINDERZ" --inner-cylinder-radius "$CYLINDERRADIUS" --inner-cylinder-length "$CYLINDERLENGHT
else
    echo "Will not use cylinder" 
    INNERCYLINDERSETTINGS=""
fi

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
echo "CYLINDER LINE    : ""\""$INNERCYLINDERSETTINGS"\""


OUTNAME=MuonGun_step1_${RUNNR}_${FILE_NR}.i3.bz2
OUTDIR=/home/dvirhilu/projects/rpp-kenclark/dvirhilu/P_ONE_dvirhilu/I3Files/generated/muongun_step1
echo "OUTFILE NAME : " ${OUTNAME}
$i3env python $script -o ${OUTDIR}/${OUTNAME} -g ${GCD_FILE} -d ${RUNNR} --f ${FILE_NR} --numevents ${NUMEVENTS} --min-energy ${MINENERGY} --max-energy ${MAXENERGY} --power-law-index ${POWERLAWINDEX} --power-law-offset ${POWERLAWOFFSET} ${INNERCYLINDERSETTINGS} ${KDESETTING}

date 