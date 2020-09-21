#ifndef SIMPLEUTIL_H_
#define SIMPLEUTIL_H_
#include <stdio.h> 
#include <stdlib.h> 

#include "parameters.h"
#include "mmio.h"

extern void print_mat( char title[], double *A, int m, int n );
extern void print_first( char title[], double *A, int m, int n, int d );
extern double* read_mat(const char * fn);
extern struct size_m get_size(const char * fn);

#endif /*SIMPLEUTIL_H_*/
