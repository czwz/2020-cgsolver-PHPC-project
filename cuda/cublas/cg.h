#ifndef CG_H_
#define CG_H_

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <math.h>

extern "C" {
    #include "util.h"
    #include "parameters.h"
    #include <cblas.h>
}

//extern void cgsolver( double *A, double *b, double *x, int m, int n );
extern "C" double * init_source_term(int n, double h);
extern "C" void cgsolver( double *A, double *b, double *x, int m, int n );
extern "C" void smvm(int m, const double* val, const int* col, const int* row, const double* x, double* y);
extern "C" void cgsolver_sparse( double *Aval, int *, int *, double *b, double *x, int n);

#endif /*CG_H_*/
