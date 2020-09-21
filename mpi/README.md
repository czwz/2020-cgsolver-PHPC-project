This folder contains the MPI implementation for cg-solver.
=====================

1.HOW TO COMPILE?

 - separately modify the c-source code, e.g. cg.c/ cg_blas.c/ etc

 - $module load gcc mvapich2 openblas
 - $mpicc -I/$OPENBLAS_ROOT/include MODIFIED-FILE.c -c -lcuda
 - $make

2.HOW TO RUN?

 - reserve nodes on cluster
 - $srun ./conjugategradient lap2D_5pt_n100.mtx

 - alternatively, job script is provided: script.sh,

```	
#!/bin/bash -l
#SBATCH --reservation cours-phpc2020
#SBATCH --account phpc2020
#SBATCH --partition parallel
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 25
#SBATCH --cpus-per-task 1

srun ./conjugategradient lap2D_5pt_n100.mtx
```

- $sbatch script
