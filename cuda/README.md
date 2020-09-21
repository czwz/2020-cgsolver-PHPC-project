This folder contains the CUDA implementation for cg-solver.
It do the parallisation based on cublas library provided by NVIDIA.

As for metrix-vector multiplication,
(i) cublas: implement of cublase dgemv
(ii) naive: code from scratch

=====================

1.HOW TO COMPILE?

 - separately modify the c-source code, e.g. cg.c/ cg_blas.c/ etc

 - $source /ssoft/spack/bin/slmodules.sh -s  x86_E5v2_Mellanox_GPU
 - $module load nvcc openblas
 - $nvcc -I/$OPENBLAS_ROOT/include MODIFIED-FILE.c -c -lcuda
 - $make

2.HOW TO RUN?

 - reserve nodes on cluster
 - $srun ./conjugategradient lap2D_5pt_n100.mtx

 - alternatively, job script is provided: script.sh,


```	
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
```

- $sbatch script

