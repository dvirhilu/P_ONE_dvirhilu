#!/bin/zsh
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --output=/home/dvirhilu/scratch/sbatchLogFiles/FluxData/arrayjob_%A_%a.log
#SBATCH --account=rpp-kenclark
#SBATCH --job-name=FluxData
date
echo "Sleeping for 30+-10 seconds to avoid hammering filesystems" 
sleep $((20 + RANDOM % 20))

startsecond=$(date +%s)
echo "Start second: " $startsecond 

# module purge
module --force purge

echo "This is step 1 of MuonGun simulations! "
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`

echo "All arguments: " $@ 
sleep 2
echo "Starting the singularity job"
singularity exec --bind /cvmfs --bind /scratch/dvirhilu --bind /project/6008051/dvirhilu --bind /project/6008051/hignight --bind /home/dvirhilu --bind /home/hignight --nv /project/6008051/hignight/singularity_images/centos7.img /project/6008051/dvirhilu/P_ONE_dvirhilu/src/I3FileCode/simAnalysis/fluxData/getFluxData_job.sh $@
date

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
