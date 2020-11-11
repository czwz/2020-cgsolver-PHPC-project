# PHPC - CONJUGATE GRADIENT PROJECT

HOWTO COMPILE AND RUN

```diff
-This version only works for regular matrix (not sparse matrix).
-This version is still under development.
```

Requirements :

- a recent compiler (ex: gcc or intel)
- a BLAS library (ex: openblas or intel MKL)
- a MPI library (ex: mvapich2)

compile on your machine (example for mpi) :

```
$ module load gcc openblas mvapich2
$ git clone https://github.com/tim010007/2020-cgsolver-PHPC-project.git
$ cd phpc_final-project_slwang/mpi
$ make
```

 run by provided script (job.sh) :
 1.serial version:
```
$ ./conjugategradient lap2D_5pt_n100.mtx
```
 2.parallel version (mpi/ cuda) on cluster
```
$ sbatch job
```
Output should be like (in slurm file):

```
cat slurm-<id>

Matrix loaded from file lap2D_5pt_n100.mtx
N = 10000
nz = 29800
val[0] = 4.000000
Call cgsolver() on matrix size (10000 x 10000)
        [STEP 488] residual = 1.103472E-10
Time for CG (dense solver)  = 63.611155 [s]


```

The given example is a 5-points stencil for the 2D Laplace problem. The matrix is in sparse format.
The matrix format is [Matrix Market format (.mtx extension)](https://sparse.tamu.edu/). You can use other matrices there or create your own.
