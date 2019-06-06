#!/bin/zsh
#SBATCH --time=23:59:00
#SBATCH --mem=8G
#SBATCH --output=/scratch/terliuk/job_outputs/muongun/step2/log_mg_s2.%A_%a.log
#SBATCH --account=rpp-dgrant
#SBATCH --job-name=MG_step2_%A
#SBATCH --gres=gpu:1
date
echo "Sleeping for 30+-10 seconds to avoid hammering filesystems" 
sleep $((20 + RANDOM % 20))

startsecond=$(date +%s)
echo "Actual start second: " $startsecond 
# source /home/terliuk/scripts/muongun_scripts/mod_purge.sh
# source /home/hignight/setup_oscnext.sh
echo "This is step 2 of MuonGun simulations! "
echo "1) opencl_vendor_path:" ${OPENCL_VENDOR_PATH}
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.1.1/setup.sh`
echo "2) opencl_vendor_path: " ${OPENCL_VENDOR_PATH}
# generating name of the output file
echo "SLURM TASK ID : " $SLURM_ARRAY_TASK_ID
echo "All arguments: " $@ 
sleep 2
echo "Starting the singularity job"
singularity exec --bind /tmp --bind /cvmfs --bind /scratch/terliuk --bind /scratch/hignight --bind /project/6008051/terliuk --bind /project/6008051/hignight --bind /home/terliuk --bind /home/hignight --nv /project/6008051/hignight/singularity_images/centos7.img /home/terliuk/scripts/muongun_scripts/MG_step2_job.sh $@

date

endsecond=$(date +%s)
echo "End second: " $endsecond 
echo "This job took : "`expr $endsecond - $startsecond`" s"
