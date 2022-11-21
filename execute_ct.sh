#!/bin/bash
#SBATCH -J build_tree        # Name of job
#SBATCH -A b1042             # Allocation, could switch to b1059 if need longer wall time
#SBATCH -p genomicsguestA    # Queue
#SBATCH -t 48:00:00          # Walltime/duration of the job
#SBATCH --cpus-per-task=4    # Number of Cores (Processors, cpus) for each task
#SBATCH --ntasks=1           # maybe set to 1 b/c we want to run one step at a time? Original = 3
#SBATCH --mem-per-cpu=25G    # Memory per node in GB needed for a job. Also see --mem-per-cpu

# set the date as parameter
radius=$1

# run it
ct *.jpg  --radius ${radius}  -d --header > results.txt
