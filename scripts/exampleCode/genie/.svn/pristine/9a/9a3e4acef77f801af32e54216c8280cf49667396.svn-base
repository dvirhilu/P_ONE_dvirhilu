#!/bin/zsh
#SBATCH --time=11:59:00
#SBATCH --mem=8G
#SBATCH --output=/home/hignight/scratch/mass_production_official/log_level12/level12_s3.%A_%a.log
#SBATCH --account=rpp-kenclark
#SBATCH --job-name=genie_level12
date
echo "Sleeping for 30+-10 seconds to avoid hammering filesystems" 
sleep $((20 + RANDOM % 20))
startsecond=$(date +%s)
echo "Start second: " $startsecond 


eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "This is step 3 of neutrino simulations! "
echo "SLURM TASK ID : " $SLURM_ARRAY_TASK_ID
echo "All arguments: " $@ 
sleep 2
echo "Starting the singularity job"
singularity exec --bind /tmp --bind /cvmfs --bind /scratch/terliuk --bind /scratch/hignight --bind /project/6008051/terliuk --bind /project/6008051/hignight --bind /project/6033576/terliuk --bind /project/6033576/hignight --bind /home/terliuk --bind /home/hignight --nv /project/6008051/hignight/singularity_images/centos7.img /project/6008051/hignight//mass_production_official/oscnext_scripts/genie_level12_jobs.sh $@
print "All done"

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
date


