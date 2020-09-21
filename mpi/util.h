#ifndef SIMPLEUTIL_H_
#define SIMPLEUTIL_H_
#include <stdio.h> 
#include <stdlib.h> 
#include "parameters.h"
#include "mmio.h"


void print_mat( char title[], double *A, int m, int n );
void print_first( char title[], double *A, int m, int n, int d );
double* read_mat(const char * restrict fn);
struct size_m get_size(const char * restrict fn);

#endif /*SIMPLEUTIL_H_*/

