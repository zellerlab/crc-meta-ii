#!/bin/bash 
# SLURM job settings
#SBATCH -A your_project_account       
#SBATCH --time=2-00:00:00             
#SBATCH --mem=64G                     
#SBATCH --cpus-per-task=1          
#SBATCH --ntasks=1               
#SBATCH --array=1-400               
#SBATCH --job-name=get_shap          

# ------------------------------
# Script for running SHAP testing in array mode
# Assumes `SHAP_job_list.sh` contains one Rscript command per line.
# Adjust the paths below as needed.
# ------------------------------

# Navigate to working directory (to the crc-meta-ii/src/analysis folder)
cd /path/to/your/shap/scripts      # <- REPLACE with actual path 

# Extract the command corresponding to this SLURM array index and execute it
cat SHAP_job_list.sh |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
