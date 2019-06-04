#!/bin/zsh
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --output=/home/dvirhilu/scratch/sbatchLogFiles/genie_step_2/arrayjob_%A_%a.log
#SBATCH --account=rpp-kenclark
#SBATCH --job-name=genie_step2_%A
#SBATCH --gres=gpu:1
date
echo "Sleeping for 30+-10 seconds to avoid hammering filesystems" 
sleep $((20 + RANDOM % 20))
startsecond=$(date +%s)
echo "Actual start second: " $startsecond 

echo "This is step 2 of neutrino simulations! "
echo "1) opencl_vendor_path:" ${OPENCL_VENDOR_PATH}
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "2) opencl_vendor_path: " ${OPENCL_VENDOR_PATH}
# generating name of the output file
echo "print =========================================="
echo "print SLURM_JOB_ID = $SLURM_JOB_ID"
echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "print =========================================="

echo "SLURM TASK ID : " $SLURM_ARRAY_TASK_ID
echo "All arguments: " $@ 
sleep 2
echo "purging modules"
module --force purge
echo "Starting the singularity job"
singularity exec --bind /tmp --bind /cvmfs --bind /scratch/dvirhilu --bind /scratch/hignight --bind /project/6008051/dvirhilu --bind /project/6008051/hignight --bind /home/dvirhilu --bind /home/hignight --nv /project/6008051/hignight/singularity_images/centos7.img /project/6008051/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/genie/genie_step2_job.sh $@

date

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
