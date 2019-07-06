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
script=/project/6008051/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/genie/step_2_clsim_setCrossE.py
echo "Will use script: " $script

RUNTYPE=$1
echo "Run type: " $RUNTYPE
FLV=$2
echo "Flavor is: " ${FLV}
E=$3
echo "Energy Range is: " ${E}
MEDIUM=$4
echo "Medium will be: " ${MEDIUM}

if [ "$RUNTYPE" == "testString" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=TestString_n15_b100.0_v50.0_l1_simple_spacing
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/testStrings/${GCDNAME}.i3.gz
 
elif [ "$RUNTYPE" == "HorizGeo" ]; then
    echo "Found configuration for " $RUNTYPE
    GCDNAME=HorizGeo_n10_b100.0_a90.0_l1_linear_reset_offset_exp_r_spacing
    echo "Name of GCD File Used: " $GCDNAME
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/uncorHorizGeo/${GCDNAME}.i3.gz
 
elif [ "$RUNTYPE" == "IceCube" ]; then
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
    echo "Using IceCube GCD"
    
else 
    echo "No configuration for " $RUNTYPE "... exiting"
    exit
fi

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

CROSS_E=30 
EFF=1.0
RUNNUM=${NU}0000
echo "Run Number        : "${RUNNUM}
echo "Flavor            : "${FLV}
echo "Energy Range      : "${E}
echo "medium model      : "${MEDIUM}
echo "cross-over E      : "${CROSS_E}
echo "DOM eff UnshadowedFraction: "${EFF}

INNAME=${FLV}_${E}_${FILE_NR}_step1.i3.zst
INDIR=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/genie/genie_step1

OUTNAME=${FLV}_${E}_${FILE_NR}_step2.i3.zst
OUTDIR=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/genie/genie_step2

echo "INNAME: " ${INDIR}/${FLV}/${INNAME}
echo "OUTNAME: " ${OUTDIR}/${FLV}/${OUTNAME}

echo "GCD: " $GCD_FILE

$i3env python $script -t -i ${INDIR}/${FLV}/${INNAME} -g ${GCD_FILE} -o ${OUTDIR}/${FLV}/${OUTNAME} -r ${RUNNUM} -l ${FILE_NR} -c ${CROSS_E} -e ${EFF} -m $MEDIUM
