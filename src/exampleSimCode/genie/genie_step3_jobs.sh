#!/bin/bash
date

echo "Starting the neutrino job"
echo "Argument line : " $@

echo "TASK ID " $SLURM_ARRAY_TASK_ID
FILE_NR=`expr $SLURM_ARRAY_TASK_ID - 1`
FILE_NR=`printf "%06d\n" $FILE_NR`
echo "Filename ID : " $FILE_NR

eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

i3env=/home/hignight/work/oscNext_official/oscNext/build_trunk_jan21_py2_v3.1.1/env-shell.sh
script=/project/6008051/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/genie/step_3_det_general_lowdt.py
I3_SRC=/home/hignight/work/oscNext_official/oscNext/trunk

echo "Will use i3 environment: " ${i3env}
echo "Will use I3_SRC : " ${I3_SRC}
echo "Will use script: " $script


RUNNUM=$1
echo "Run number: " $RUNNUM
FLV=$2
echo "Flavor is: " ${FLV}
E=$3
echo "Energy Range is: " ${E}


case ${RUNNUM} in
    0000)
	#all these are default values
	ICEMODEL=spice_3.2.1
	CROSS_E=30 
	DOMEFF=1.0
	HOLEICE=as.flasher_p1_0.30_p2_-1
	;;
    *)
	echo "No configuration for " $RUNNUM "... exiting"
	exit 1
	;;
esac

#Get FLV setup
case ${FLV} in
    NuE)
        NU=12
        ;;
    NuMu)
        NU=14
        ;;
    NuTau)
        NU=16
        ;;
    *)
        echo ${FLV} " is not an acceptable neutrino type (NuE NuMu NuTau)"
        exit 3
        ;;

esac

RUNNUM=${NU}${RUNNUM}

INNAME=${FLV}_${E}_${RUNNUM}_${FILE_NR}_step2.i3.zst
INDIR=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/genie_step2


OUTNAME=${FLV}_${E}_${RUNNUM}_${FILE_NR}_step3.i3.zst
OUTDIR=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/generated/genie_step3

echo "INNAME: " ${INDIR}/${FLV}/${INNAME}
echo "OUTNAME: " ${OUTDIR}/${FLV}/${OUTNAME}

echo "ICE MODEL      : " $ICEMODEL
echo "DOM EFFICINECY : " $DOMEFF
echo "HOLE ICE       : " $HOLEICE

GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
echo "GCD: " $GCD_FILE

echo "Starting the job"
$i3env python $script -i ${INDIR}/$INNAME -g $GCD_FILE -o ${OUTDIR}/${OUTNAME} -r ${RUNNUM} -f ${FILE_NR} -e ${DOMEFF}  -l ${HOLEICE} -m ${ICEMODEL}

