#!/bin/bash -l
#SBATCH --reservation cours-phpc2020
#SBATCH --account phpc2020
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 20
#SBATCH --cpus-per-task 1

module purge
module load intel
module load intel-mpi
srun conjugategradient lap2D_5pt_n100.mtx
#perf record -o perf.data --call-graph dwarf srun conjugategradient lap2D_5pt_n100.mtx
