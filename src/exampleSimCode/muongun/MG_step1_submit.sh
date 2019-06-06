#!/bin/zsh
#SBATCH --time=05:59:00
#SBATCH --mem=4G
#SBATCH --output=/scratch/terliuk/job_outputs/muongun/step1/log_mg_s1.%A_%a.log
#SBATCH --account=rpp-dgrant
#SBATCH --job-name=MG_step1
date
echo "Sleeping for 30+-10 seconds to avoid hammering filesystems" 
sleep $((20 + RANDOM % 20))

startsecond=$(date +%s)
echo "Start second: " $startsecond 

# source /home/terliuk/scripts/muongun_scripts/mod_purge.sh
# source /home/hignight/setup_oscnext.sh
echo "This is step 1 of MuonGun simulations! "
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
singularity exec --bind /cvmfs --bind /scratch/terliuk --bind /scratch/hignight --bind /project/6008051/terliuk --bind /project/6008051/hignight --bind /home/terliuk --bind /home/hignight --nv /project/6008051/hignight/singularity_images/centos7.img /home/terliuk/scripts/muongun_scripts/MG_step1_job.sh $@
date

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
