#!/bin/bash

# Job name:
#SBATCH --job-name=run_match
#
# Project:
#SBATCH --account=ec332
#
# Wall time limit:
#SBATCH --time=00-04:00:00

# Resources
#SBATCH --mem=84G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32

script_dir="/fp/homes01/u01/ec-karinal/New_bubbles"
script_name="event_matching_main.py"

yr_mn=$1
echo $yr_mn

## Set up job environment:
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

# Set the ${PS1} (needed in the source of the Anaconda environment)
export PS1=\$

module --quiet purge  # Reset the modules to the system default
module load Miniconda3/22.11.1-1
module list

# Set the python environment
source ${EBROOTMINICONDA3}/etc/profile.d/conda.sh
conda deactivate
conda activate /fp/homes01/u01/ec-karinal/.conda/envs/saff

## Do some work:
srun python ${script_dir}/${script_name} "$yr_mn"
