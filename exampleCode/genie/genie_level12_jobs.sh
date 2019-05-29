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

INDIR=/project/6008051/hignight/mass_production_official/step3/${FLV}
OUTDIR=/project/6008051/hignight/mass_production_official/level12/${FLV}

OUTNAME1=${FLV}_${RUNNUM}_${FILE_NR}_level1.zst
OUTNAME2=${FLV}_${RUNNUM}_${FILE_NR}_level2.zst

echo "OUTNAME: " ${OUTDIR}/${OUTNAME1}

echo "ICE MODEL      : " $ICEMODEL
echo "DOM EFFICINECY : " $DOMEFF
echo "HOLE ICE       : " $HOLEICE

GCD_FILE=/project/6008051/hignight/GCD_with_noise/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_withScaledNoise.i3.gz
echo "GCD: " $GCD_FILE

echo "Starting the job"
echo "--- Level 1---"
level1script=$I3_SRC/filterscripts/resources/scripts/SimulationFiltering.py
echo "Level1 script : " $level1script
echo $i3env $level1script -g $GCD_FILE -i ${INDIR}/${FLV}_A_${RUNNUM}_${FILE_NR}_step3.zst,${INDIR}/${FLV}_B_${RUNNUM}_${FILE_NR}_step3.zst,${INDIR}/${FLV}_C_${RUNNUM}_${FILE_NR}_step3.zst,${INDIR}/${FLV}_D_${RUNNUM}_${FILE_NR}_step3.zst -o ${OUTDIR}/${OUTNAME1} --photonicsdir=/cvmfs/icecube.opensciencegrid.org/data/photon-tables/

$i3env $level1script -g $GCD_FILE -i ${INDIR}/${FLV}_A_${RUNNUM}_${FILE_NR}_step3.zst,${INDIR}/${FLV}_B_${RUNNUM}_${FILE_NR}_step3.zst,${INDIR}/${FLV}_C_${RUNNUM}_${FILE_NR}_step3.zst,${INDIR}/${FLV}_D_${RUNNUM}_${FILE_NR}_step3.zst -o ${OUTDIR}/${OUTNAME1} --photonicsdir=/cvmfs/icecube.opensciencegrid.org/data/photon-tables/
echo "--- Level 2---"
level2script=$I3_SRC/filterscripts/resources/scripts/offlineL2/process.py
echo "Level2 script : " $level2script
echo $i3env $level2script -s -g $GCD_FILE -i ${OUTDIR}/${OUTNAME1} -o ${OUTDIR}/${OUTNAME2} --photonicsdir=/cvmfs/icecube.opensciencegrid.org/data/photon-tables/
$i3env $level2script -s -g $GCD_FILE -i ${OUTDIR}/${OUTNAME1} -o ${OUTDIR}/${OUTNAME2} --photonicsdir=/cvmfs/icecube.opensciencegrid.org/data/photon-tables/

rm ${OUTDIR}/${OUTNAME1}
