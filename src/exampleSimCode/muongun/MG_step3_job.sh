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
I3_SRC=/home/hignight/work/oscNext_official/oscNext/trunk
echo "Will use i3 environment: " ${i3env}
echo "Will use I3_SRC : " ${I3_SRC}
script=/projects/6008051/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/muongun/step_3_det_general_lowdt.py
echo "Will use script: " $script

OUTRUN=$1

if [ "$OUTRUN" == "139005" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139005
    ICEMODEL="spice_3.2.1"
    DOMEFF=1.00
    HOLEICE="as.flasher_p1_0.30_p2_-1"
elif [ "$OUTRUN" == "139006" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139006
    ICEMODEL="spice_3.2.1"
    DOMEFF=1.00
    HOLEICE="as.flasher_p1_0.30_p2_-1"
elif [ "$OUTRUN" == "139008" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139008
    ICEMODEL="spice_3.2.1"
    DOMEFF=1.00
    HOLEICE="as.flasher_p1_0.30_p2_-1"
elif [ "$OUTRUN" == "139010" ]; then
    echo "Found configuration for " $OUTRUN
    INRUN=139010
    ICEMODEL="spice_3.2.1"
    DOMEFF=1.00
    HOLEICE="as.flasher_p1_0.30_p2_-1"
else 
    echo "No configuration for " $OUTRUN "... exiting"
    exit
fi

echo "IN RUN         : " $INRUN
echo "OUT RUN        : " $OUTRUN
echo "ICE MODEL      : " $ICEMODEL
echo "DOM EFFICINECY : " $DOMEFF
echo "HOLE ICE       : " $HOLEICE

INFILENAME=MuonGun_step2_${INRUN}_${FILE_NR}.i3.bz2
INFOLDER=/projects/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/muongun_step2
echo "INFILEPATH: " ${INFOLDER}/$INFILENAME
OUTFILENAME=MuonGun_step3_${OUTRUN}_${FILE_NR}.i3.bz2
echo "OUTFILENAME : " $OUTFILENAME
GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
echo "GCD: " $GCD_FILE

ls -ltrk /scratch/terliuk/sim/muongun/
echo "Starting the job"
$i3env python $script -i ${INFOLDER}/$INFILENAME -g $GCD_FILE -o /scratch/terliuk/sim/muongun/tmpStep3_${OUTFILENAME} -r $OUTRUN -f $SLURM_ARRAY_TASK_ID -e $DOMEFF  -l $HOLEICE -m $ICEMODEL
ls -ltrk /scratch/terliuk/sim/muongun/
echo "-----*** MOVING to Level1 and Level2 ***------"
echo "--- Level 1---"
level1script=$I3_SRC/filterscripts/resources/scripts/SimulationFiltering.py
echo "Level1 script : " $level1script
echo "Level1 output : " /scratch/terliuk/sim/muongun/tmpLevel1_${OUTFILENAME}
$i3env $level1script -g $GCD_FILE -i /scratch/terliuk/sim/muongun/tmpStep3_$OUTFILENAME -o /scratch/terliuk/sim/muongun/tmpLevel1_${OUTFILENAME} --photonicsdir=/cvmfs/icecube.opensciencegrid.org/data/photon-tables/
ls -ltrk /scratch/terliuk/sim/muongun/
echo "--- Level 2---"
level2script=$I3_SRC/filterscripts/resources/scripts/offlineL2/process.py
echo "Level2 script : " $level2script
$i3env $level2script -s -g $GCD_FILE -i /scratch/terliuk/sim/muongun/tmpLevel1_$OUTFILENAME -o /scratch/terliuk/sim/muongun/$OUTFILENAME --photonicsdir=/cvmfs/icecube.opensciencegrid.org/data/photon-tables/
ls -ltrk /scratch/terliuk/sim/muongun/
mv /scratch/terliuk/sim/muongun/$OUTFILENAME /home/terliuk/projects/rpp-dgrant/terliuk/simulations/muongun/Level2/$OUTRUN/$OUTFILENAME 
ls -ltrk /scratch/terliuk/sim/muongun/
rm /scratch/terliuk/sim/muongun/tmpStep3_$OUTFILENAME /scratch/terliuk/sim/muongun/tmpLevel1_$OUTFILENAME
