#include "util.h"




// ====================================================================
// Implementation of useful functions
// ====================================================================

/*
print a matrix A of size m x n
*/

void print_mat( char title[], double *A, int m, int n )
{
	double x;
	printf("%s",title);

	for ( int i = 0; i < m; i++ ){
		for ( int j = 0; j < n; j++ ){
			x = A[id(m,i,j)];   //if ( abs( x ) < NEARZERO ) x = 0.0;
			printf("%f (%d)\t",x,id(m,i,j));
		}
		printf("\n");
	}
}

/*
print first d elements from a matrix A of size m x n
*/

void print_first( char title[], double *A, int m, int n, int d )
{
	double x;
	printf("%s",title);

	for ( int i = 0; i < d; i++ ){
		x = A[i];   //if ( abs( x ) < NEARZERO ) x = 0.0;
		printf("%f (%d)\t",x,i);
	}
	printf("\n");
}



/*
reads a matrix A of size m x n

inspired by https://math.nist.gov/MatrixMarket/mmio/c/example_read.c
*/

double* read_mat(const char * restrict fn)
{
	double * mat;
	int n,m;

	int ret_code;
	MM_typecode matcode;
	FILE *f;
	int nz;   
	int i, *I, *J;
	double tmp;

	if ((f = fopen(fn, "r")) == NULL) 
	{
		printf("Could not open matrix");
		exit(1);
	}

	if (mm_read_banner(f, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}


//	Matrix is sparse
	if (mm_is_matrix(matcode) && mm_is_coordinate(matcode) ) {
		if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nz)) !=0)
		{
			exit(1);
		}else{
			if ( (m==1) || (n==1) ){
				printf("size of vector = %d x %d\n",m,n);
			}else{
				printf("size of matrix = %d x %d\n",m,n);
			}
		}
	}
//	Matrix is not sparse
	else if (mm_is_matrix(matcode) && mm_is_array(matcode) )
	{
//		printf("case : %%MatrixMarket matrix array real symmetric \n");
		if ((ret_code = mm_read_mtx_array_size(f, &m, &n)) !=0)
		{
			exit(1);
		}else{
			if ( (m==1) || (n==1) ){
				printf("size of vector = %d x %d\n",m,n);
			}else{
				printf("size of matrix = %d x %d\n",m,n);
			}
		}
		nz = m*n;
	}
//	Matrix is something else
	else{
		printf("Sorry, this application does not support ");
		printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
		exit(1);
	}

/* reserve memory for matrices */

	I = (int *) malloc(nz * sizeof(int));
	J = (int *) malloc(nz * sizeof(int));
	mat = (double*) malloc(m*n * sizeof(double));

	for (i=0; i<m*n; i++)
	{
		mat[i] = 0.;
	}


    /*  NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */


// Matrix is symmetric

	if (mm_is_symmetric(matcode)){

		for (i=0; i<nz; i++)
		{
			if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &tmp)==3){
				I[i]--;  /* adjust from 1-based to 0-based */
				J[i]--;
				mat[I[i]*m+J[i]] = tmp;
				mat[J[i]*m+I[i]] = tmp;
			}else if(fscanf(f, "%lg\n",&tmp)==1){
				mat[i] = tmp;
			}else{
				break;
			}
		}
	}

// Matrix is not symmetric

	if (mm_is_general(matcode)){
		for (i=0; i<nz; i++)
		{
			if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &tmp)==3){
				I[i]--;  /* adjust from 1-based to 0-based */
				J[i]--;
				mat[I[i]*m+J[i]] = tmp;
			}else if(fscanf(f, "%lg\n",&tmp)==1){
				mat[i] = tmp;
			}else{
				break;
			}
		}
	}


	if (f !=stdin) fclose(f);

//	print_mat( "matrice : \n", mat, m, n );

	return mat;
}



/*
Reads the size of the matrix
*/

/*
reads a matrix A of size m x n
*/

struct size_m get_size(const char * restrict fn)
{
	int n,m,nz;
	FILE *f;
	struct size_m size;
	int ret_code;
	MM_typecode matcode;

	if ((f = fopen(fn, "r")) == NULL) 
	{
		printf("Could not open matrix");
		exit(1);
	}

	if (mm_read_banner(f, &matcode) != 0)
	{
		printf("Could not process Matrix Market banner.\n");
		exit(1);
	}

	if (mm_is_matrix(matcode) && mm_is_coordinate(matcode) ) {
		if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nz)) !=0)
		{
			exit(1);
		}
	}
	else if (mm_is_matrix(matcode) && mm_is_array(matcode) )
	{
		if ((ret_code = mm_read_mtx_array_size(f, &m, &n)) !=0)
		{
			exit(1);
		}
	}


//	if ((ret_code = mm_read_mtx_crd_size(f, &m, &n, &nz)) !=0)
//		exit(1);
	size.m = m;
	size.n = n;

	if (f !=stdin) fclose(f);

	return size;
}


