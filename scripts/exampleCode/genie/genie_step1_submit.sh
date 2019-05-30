#!/bin/zsh
#SBATCH --time=11:59:00
#SBATCH --mem=4G
#SBATCH --output=/home/dvirhilu/scratch/sbatchLogFiles/genie_step_1/arrayjob_%A_%a.log
#SBATCH --account=rpp-kenclark
#SBATCH --job-name=genie_step1
date
echo "Sleeping for 30+-10 seconds to avoid hammering filesystems" 
sleep $((20 + RANDOM % 20))

startsecond=$(date +%s)
echo "Start second: " $startsecond 

echo "This is step 1 of neutrino simulations! "
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
# generating name of the output file
echo "print =========================================="
echo "print SLURM_JOB_ID = $SLURM_JOB_ID"
echo "print SLURM_JOB_NODELIST = $SLURM_JOB_NODELIST"
echo "print =========================================="

echo "SLURM TASK ID : " $SLURM_ARRAY_TASK_ID
echo "All arguments: " $@ 
sleep 2
echo "Starting the singularity job"
singularity exec --bind /cvmfs --bind /scratch/dvirhilu --bind /scratch/hignight --bind /project/6008051/dvirhilu --bind /project/6008051/hignight --bind /home/dvirhilu --bind /home/hignight --nv /project/6008051/hignight/singularity_images/centos7.img /project/6008051/dvirhilu/P_ONE_dvirhilu/scripts/exampleCode/genie/genie_step1_job.sh $@
date

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
