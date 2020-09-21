#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=1:0:0
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=gpu_free
##SBATCH --reservation cours-phpc2020
#SBATCH --account phpc2020

source /ssoft/spack/bin/slmodules.sh -s  x86_E5v2_Mellanox_GPU

module load gcc cuda
module load openblas
srun ./conjugategradient lap2D_5pt_n100.mtx

