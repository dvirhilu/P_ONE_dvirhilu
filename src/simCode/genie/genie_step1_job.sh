#!/bin/bash
date

echo "Starting the neutrino job"
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
script=/project/6008051/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/genie/step_1_genie.py
echo "Will use script: " $script

RUNTYPE=$1
echo "Run type: " $RUNTYPE
FLV=$2
echo "Flavor is: " ${FLV}
E=$3
echo "Energy Range is: " ${E}

OUTDIR=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/genie/genie_step1

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

elif [ "$RUNTYPE" == "cube" ]; then
    echo "Found configuration for " $RUNTYPE
    GCD_FILE=/project/6008051/dvirhilu/P_ONE_dvirhilu/I3Files/gcd/cube/cubeGeometry_1600_15_50.i3.gz

else 
    echo "No configuration for " $RUNTYPE "... exiting"
    exit
fi

#Get FLV setup
case ${FLV} in
    NuE)
        NU=12
        case ${E} in
            A)
                NEVENTS=4500    # modified from 450000
                ;;
            B)
                NEVENTS=1000    # modified from 100000
                ;;
            C)
                NEVENTS=1000    # modified from 100000
                ;;
            D)
                NEVENTS=575     # modified from 57500
                ;;
	    *)
		echo ${E} " is not an acceptable neutrino energy range (A B C D)"
		exit 2
		;;
        esac
        ;;
    NuMu)
        NU=14
        case ${E} in
            A)
                NEVENTS=4080    # modified from 408000
                ;;
            B)
                NEVENTS=4400    # modified from 440000
                ;;
            C)
                NEVENTS=1000     # modified from 57500
                ;;
            D)
                NEVENTS=67      # modified from 6700
                ;;
	    *)
		echo ${E} " is not an acceptable neutrino energy range (A B C D)"
		exit 2
		;;
        esac
        ;;
    NuTau)
        NU=16
        case ${E} in
            A)
                NEVENTS=3000    # modified from 300000
                ;;
            B)
                NEVENTS=3750    # modified from 375000
                ;;
            C)
                NEVENTS=5000    # modified from 200000
                ;;
            D)
                NEVENTS=260     # modified from 26000
                ;;
	    *)
		echo ${E} " is not an acceptable neutrino energy range (A B C D)"
		exit 2
		;;
        esac
        ;;
    *)
        echo ${FLV} " is not an acceptable neutrino type (NuE NuMu NuTau)"
        exit 3
        ;;

esac

RUNNUM=${NU}0000
echo "Run Number        : "${RUNNUM}
echo "Flavor            : "${FLV}
echo "Energy Range      : "${E}
echo "Number of Events  : "${NEVENTS}


OUTNAME=${FLV}_${E}_${FILE_NR}_step1.i3.zst

echo "OUTFILE NAME : " ${OUTNAME}
$i3env python $script -o ${OUTDIR}/${FLV}/${OUTNAME} -l ${FILE_NR}  -r ${RUNNUM} -n ${NEVENTS} -f ${FLV} --energy-range ${E} 

date


###########################
####     Exit codes    ####
#### 1 - No set number ####
#### 2 - No E code     ####
#### 3 - No FLV code   ####
###########################

