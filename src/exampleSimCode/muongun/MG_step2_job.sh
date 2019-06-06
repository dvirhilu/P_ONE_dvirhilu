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
script=/home/terliuk/projects/rpp-dgrant/terliuk/LE_simulation_scripts/step_2_clsim_setCrossE.py
echo "Will use script: " $script

OUTRUN=$1

if [ "$OUTRUN" == "139005" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139005
    ICEMODEL="spice_3.2.1"
    CROSSENERGY=200.0
elif [ "$OUTRUN" == "139006" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139006
    ICEMODEL="spice_3.2.1"
    CROSSENERGY=200.0
elif [ "$OUTRUN" == "139008" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139008
    ICEMODEL="spice_3.2.1"
    CROSSENERGY=200.0
elif [ "$OUTRUN" == "139010" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139010
    ICEMODEL="spice_3.2.1"
    CROSSENERGY=200.0
else 
    echo "No configuration for " $OUTRUN "... exiting"
    exit
fi
echo "INRUN       : " $INRUN
echo "ICEMODEL    : " $ICEMODEL
echo "CROSSENERGY : " $CROSSENERGY
INFILENAME=MuonGun_step1_${INRUN}.${FILE_NR}.i3.bz2
INFOLDER=/project/6008051/terliuk/simulations/muongun/step1/$INRUN/

echo "INFILEPATH: " ${INFOLDER}/$INFILENAME
OUTFILENAME=MuonGun_step2_${OUTRUN}.${FILE_NR}.i3.bz2
echo "OUTFILENAME : " $OUTFILENAME
GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
echo "GCD: " $GCD_FILE

ls -ltrk /scratch/terliuk/sim/muongun/
$i3env python $script -t -i ${INFOLDER}/$INFILENAME -g $GCD_FILE -o /scratch/terliuk/sim/muongun/$OUTFILENAME -r ${OUTRUN} -l $SLURM_ARRAY_TASK_ID -m $ICEMODEL -c $CROSSENERGY
ls -ltrk /scratch/terliuk/sim/muongun/
mv /scratch/terliuk/sim/muongun/$OUTFILENAME /project/6008051/terliuk/simulations/muongun/step2/$OUTRUN/$OUTFILENAME 
