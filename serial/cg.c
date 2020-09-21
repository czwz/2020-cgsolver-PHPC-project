#include "cg.h"

const double TOLERANCE = 1.0e-10;

/*
	cgsolver solves the linear equation A*x = b where A is 
	of size m x n

Code based on MATLAB code (from wikipedia ;-)  ):

function x = conjgrad(A, b, x)
    r = b - A * x;
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = A * p;
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        rsold = rsnew;
    end
end

*/
void cgsolver( double *A, double *b, double *x, int m, int n ){
	double * r;
	double * rt;
	double * p;
	double * pt;
	double rsold;
	double rsnew;
	double * Ap;
	double * tmp;
	double alpha;

	int incx = 1;
	int incy = 1;
	int lda = m;
	double al = 1.;
	double be = 0.;

	int k = 0;

	r = (double*) malloc(n* sizeof(double));
	rt = (double*) malloc(n* sizeof(double));
	p = (double*) malloc(n* sizeof(double));
	pt = (double*) malloc(n* sizeof(double));
	Ap = (double*) malloc(n* sizeof(double));
	tmp = (double*) malloc(n* sizeof(double));

//    r = b - A * x;
	cblas_dgemv (CblasColMajor, CblasNoTrans, m, n, al , A, lda, x, incx,  be, Ap, incy);
	cblas_dcopy(n,b,1,tmp,1);	
	cblas_daxpy(n, -1.  , Ap,1, tmp, 1 );
	cblas_dcopy(n,tmp,1,r,1);

//    p = r;
	cblas_dcopy (n, r, incx, p, incy);

	cblas_dcopy(n,r,1,rt,1);	
	cblas_dcopy(n,p,1,pt,1);	
//    rsold = r' * r;
	rsold = cblas_ddot (n,r,incx,rt,incy);

//    for i = 1:length(b)
	while ( k < n  ){
//        Ap = A * p;
		cblas_dgemv (CblasColMajor, CblasNoTrans, m, n, al, A, lda, p, incx, be, Ap, incy);
//        alpha = rsold / (p' * Ap);
		alpha = rsold / fmax( cblas_ddot(n, pt, incx, Ap, incy ), NEARZERO );
//        x = x + alpha * p;
		cblas_daxpy(n, alpha, p, 1, x, 1);
//        r = r - alpha * Ap;
		cblas_daxpy(n, -alpha, Ap, 1, r, 1);
//        rsnew = r' * r;
		rsnew = cblas_ddot (n,r,incx,r,incy);
//        if sqrt(rsnew) < 1e-10
//              break;
		if ( sqrt(rsnew) < TOLERANCE ) break;             // Convergence test
//        p = r + (rsnew / rsold) * p;
		cblas_dcopy(n,r,1,tmp,1);	
		cblas_daxpy(n, (double)(rsnew/rsold), p, 1, tmp, 1);
		cblas_dcopy(n,tmp,1,p,1);
	
		cblas_dcopy(n,p,1,pt,1);	
//        rsold = rsnew;
		rsold = rsnew;
		k++;
	}

	printf("\t[STEP %d] residual = %E\n",k,sqrt(rsold));

	free(r);
	free(rt);
	free(p);
	free(pt);
	free(Ap);
	free(tmp);
}

/*
Sparse version of the cg solver
*/

void cgsolver_sparse( double *Aval, int *Irn, int *Jcn, double *b, double *x, int n){
	double * r;
	double * rt;
	double * p;
	double * pt;
	double rsold;
	double rsnew;
	double * Ap;
	double * tmp;
	double alpha;

	int incx = 1;
	int incy = 1;

	int k = 0;

	r = (double*) malloc(n* sizeof(double));
	rt = (double*) malloc(n* sizeof(double));
	p = (double*) malloc(n* sizeof(double));
	pt = (double*) malloc(n* sizeof(double));
	Ap = (double*) malloc(n* sizeof(double));
	tmp = (double*) malloc(n* sizeof(double));

//    r = b - A * x;
	smvm(n, Aval, Jcn, Irn, x, Ap);
	cblas_dcopy(n,b,1,tmp,1);	
	cblas_daxpy(n, -1.  , Ap,1, tmp, 1 );
	cblas_dcopy(n,tmp,1,r,1);

//    p = r;
	cblas_dcopy (n, r, incx, p, incy);

	cblas_dcopy(n,r,1,rt,1);	
	cblas_dcopy(n,p,1,pt,1);	
//    rsold = r' * r;
	rsold = cblas_ddot (n,r,incx,rt,incy);



//    for i = 1:length(b)
	while ( k < n ){
//        Ap = A * p;
		smvm(n, Aval, Jcn, Irn, p, Ap);
//        alpha = rsold / (p' * Ap);
		alpha = rsold / fmax( cblas_ddot(n, pt, incx, Ap, incy ), NEARZERO );
//        x = x + alpha * p;
		memset(x, 0, n*sizeof(double));
		cblas_daxpy(n, alpha, p, 1, x, 1);
//        r = r - alpha * Ap;
		cblas_daxpy(n, -alpha, Ap, 1, r, 1);
//        rsnew = r' * r;
		rsnew = cblas_ddot (n,r,incx,r,incy);
//        if sqrt(rsnew) < 1e-10
//              break;
		if ( sqrt(rsnew) < TOLERANCE ) break;             // Convergence test
//        p = r + (rsnew / rsold) * p;
		memset(p, 0, n*sizeof(double));
		cblas_dcopy(n,r,1,tmp,1);	
		cblas_daxpy(n, (double)(rsnew/rsold), p, 1, tmp, 1);
		cblas_dcopy(n,tmp,1,p,1);
	
		cblas_dcopy(n,p,1,pt,1);	
//        rsold = rsnew;
		rsold = rsnew;
		k++;
	}

	printf("\t[STEP %d] residual = %E\n",k,sqrt(rsold));

	free(r);
	free(rt);
	free(p);
	free(pt);
	free(Ap);
	free(tmp);
}


/*
Sparse matrix vector multiplication
*/

void smvm(int m, const double* val, const int* col, const int* row, const double* x, double* y)
{
	for (int i=0; i<m; ++i) {
		y[i] = 0.0;
		for (int j=row[i]; j<row[i+1]; ++j){
			y[i] += val[j-1]*x[col[j-1]-1];
		}
	}
}




/*
Initialization of the source term b 
*/

double * init_source_term(int n, double h){
	double * f;
	int i;
	f  = (double*) malloc(n*sizeof(double*));

	for(i = 0; i < n; i++) {
		f[i] = (double)i * -2. * M_PI * M_PI * sin(10.*M_PI*i*h) * sin(10.*M_PI*i*h);
	}
	return f;
}


