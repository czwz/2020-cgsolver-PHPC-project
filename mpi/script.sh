#!/bin/bash -l
#SBATCH --reservation cours-phpc2020
#SBATCH --account phpc2020
#SBATCH --partition parallel
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 25
#SBATCH --cpus-per-task 1

srun ./conjugategradient lap2D_5pt_n100.mtx
