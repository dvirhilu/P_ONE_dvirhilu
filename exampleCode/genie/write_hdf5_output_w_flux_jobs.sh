#!/bin/bash
date

echo "Starting the neutrino job"
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

RUNNUM=$1
echo "Run number: " $RUNNUM
FLV=$2
echo "Flavor is: " ${FLV}

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

INDIR=/project/6033576/hignight/mass_production_official/level12/${FLV}
OUTDIR=/project/6033576/hignight/mass_production_official/level12/hdf5

OUTNAME=${FLV}_${RUNNUM}_${FILE_NR}_level2.hdf5

echo "OUTNAME: " ${OUTDIR}/${OUTNAME}

echo "ICE MODEL      : " $ICEMODEL
echo "DOM EFFICINECY : " $DOMEFF
echo "HOLE ICE       : " $HOLEICE

echo "Starting the job"

script=/project/6033576/hignight/mass_production_official/scripts/write_output_w_flux.py

$i3env python $script -i ${INDIR}/${FLV}_${RUNNUM}_${FILE_NR}_level2.zst -o ${OUTDIR}/${OUTNAME} 
