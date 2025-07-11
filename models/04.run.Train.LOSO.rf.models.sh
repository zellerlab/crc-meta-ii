#!/bin/bash
# SLURM job settings
#SBATCH -A your_project_account          # Replace with your SLURM project/account
#SBATCH --time=0-08:00:00                # Max wall time (adjust as needed)
#SBATCH --mem=8G                         # Memory per job
#SBATCH --cpus-per-task=10               # Number of CPUs per task
#SBATCH --ntasks=1                       # One task per job
#SBATCH --array=1-27                     # Number of jobs = number of lines in command list
#SBATCH --job-name=loso_training         # Job name

# Record the start time
start_time=$(date +%s)

# Run the command for the current array task ID
cat input.LOSO.rf.models.sh | head -n $SLURM_ARRAY_TASK_ID | tail -1 | bash

# Record the end time
end_time=$(date +%s)

# Calculate the duration in seconds
duration=$((end_time - start_time))

# Convert duration to hours, minutes, and seconds
hours=$((duration / 3600))
minutes=$(( (duration % 3600) / 60 ))
seconds=$((duration % 60))

# Print the duration
echo "Job $SLURM_ARRAY_TASK_ID completed in ${hours}h ${minutes}m ${seconds}s"

