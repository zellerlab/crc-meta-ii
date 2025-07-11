#!/bin/bash 
# SLURM job settings
#SBATCH -A your_project_account       # Replace with your SLURM project/account name
#SBATCH --time=2-00:00:00             # Max wall time (2 days)
#SBATCH --mem=64G                     # Memory per job (adjust if needed)
#SBATCH --cpus-per-task=1            # Number of CPUs per task
#SBATCH --ntasks=1                   # One task per job
#SBATCH --array=1-300               # Job array for parallel SHAP runs
#SBATCH --job-name=get_shap          # Job name

# ------------------------------
# Script for running SHAP testing in array mode
# Assumes `SHAP_job_list.sh` contains one Rscript command per line.
# Adjust the paths below as needed.
# ------------------------------

# Navigate to working directory (to the crc-meta-ii/src/analysis folder)
cd /path/to/your/shap/scripts      # <- REPLACE with actual path 

# Extract the command corresponding to this SLURM array index and execute it
cat SHAP_job_list.sh |  head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash
