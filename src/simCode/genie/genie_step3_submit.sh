#!/bin/zsh
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=/home/dvirhilu/scratch/sbatchLogFiles/genie_step_3/step3_s3.%A_%a.log
#SBATCH --account=rpp-kenclark
#SBATCH --job-name=genie_step3
date
echo "Sleeping for 30+-10 seconds to avoid hammering filesystems" 
sleep $((20 + RANDOM % 20))
startsecond=$(date +%s)
echo "Start second: " $startsecond 


eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
# generating name of the output file
echo "This is step 3 of neutrino simulations! "
echo "SLURM TASK ID : " $SLURM_ARRAY_TASK_ID
echo "All arguments: " $@ 
sleep 2
echo "purging modules"
module --force purge
echo "Starting the singularity job"
singularity exec --bind /tmp --bind /cvmfs --bind /scratch/dvirhilu --bind /scratch/hignight --bind /project/6008051/dvirhilu --bind /project/6008051/hignight --bind /home/dvirhilu --bind /home/hignight --nv /project/6008051/hignight/singularity_images/centos7.img /project/6008051/dvirhilu/P_ONE_dvirhilu/src/exampleSimCode/genie/genie_step3_jobs.sh $@
print "All done"

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
date


