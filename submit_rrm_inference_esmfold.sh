#!/bin/bash

#SBATCH --job-name=rrm_inference                     # sets the job name
#SBATCH --output=rrm_inference.%j                    # indicates a file to redirect STDOUT to; %j is the jobid. If set, must be set to a file instead of a directory or else submission will fail.
#SBATCH --error=rrm_inference.%j                     # indicates a file to redirect STDERR to; %j is the jobid. If set, must be set to a file instead of a directory or else submission will fail.
#SBATCH --time=04:00:00                                      # how long you would like your job to run; format=hh:mm:ss

#SBATCH --partition=vulcan-scavenger
#SBATCH --qos=vulcan-scavenger                                # set QOS, this will determine what resources can be requested
#SBATCH --account=vulcan-abhinav
#SBATCH --gres=gpu:rtxa5000:1

#SBATCH --nodes=1                                               # number of nodes to allocate for your job
#SBATCH --ntasks=1                                              
#SBATCH --ntasks-per-node=1                                      
#SBATCH --mem=64gb                                               # (cpu) memory required by job; if unit is not specified MB will be assumed

module load cuda
source ~/.bashrc
micromamba activate protein-lm

srun python scripts/esm_inference.py --plm_model facebook/esmfold_v1 \
                                     --msa_file data/rrm/PF00076_10000_msa_trimmed.faa \
                                     --batch_size 16 \
                                     --output_dir results/rrm/ \
                                     --job_name rrm_esmfold


wait                                                            # wait for any background processes to complete

# once the end of the batch script is reached your job allocation will be revoked