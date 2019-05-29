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
echo "Will use i3 environment: " ${i3env}
script=/project/6008051/hignight/mass_production_official/oscnext_scripts/step_2_clsim_setCrossE.py
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
	ICE=spice_3.2.1
	CROSS_E=30 
	EFF=1.2
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
echo "Run Number        : "${RUNNUM}
echo "Flavor            : "${FLV}
echo "Energy Range      : "${E}
echo "Ice model         : "${ICE}
echo "cross-over E      : "${CROSS_E}
echo "DOM eff UnshadowedFraction: "${EFF}

INNAME=${FLV}_${E}_${RUNNUM}_${FILE_NR}_step1.zst
INDIR=/project/6008051/hignight/mass_production_official/step1/

OUTNAME=${FLV}_${E}_${RUNNUM}_${FILE_NR}_step2.zst
OUTDIR=/project/6008051/hignight/mass_production_official/step2/

echo "INNAME: " ${INDIR}/${FLV}/${INNAME}
echo "OUTNAME: " ${OUTDIR}/${FLV}/${OUTNAME}

GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
echo "GCD: " $GCD_FILE

$i3env python $script -t -i ${INDIR}/${FLV}/${INNAME} -g ${GCD_FILE} -o ${OUTDIR}/${FLV}/${OUTNAME} -r ${RUNNUM} -l ${FILE_NR} -m ${ICE} -c ${CROSS_E} -e ${EFF}
