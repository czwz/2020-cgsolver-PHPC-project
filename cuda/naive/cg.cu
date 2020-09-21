#include<cublas_v2.h>

extern "C" {
    #include "cg.h"
}

const double TOLERANCE = 1.0e-10;

__device__ double atomic_Add(double* address, double val){
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                        __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}

__global__ void kernel_mv(const double * const A, const double * const X, double * const AA, int m) {
    int i,j,k,l;
    int t = blockDim.x*gridDim.x;
    int id = threadIdx.x + blockIdx.x*blockDim.x;
    double temp;

    for (i=0;i<(m*m/t);i++){
        j=(id+i*t);
        k=(id+i*t)/m;
	l=(id+i*t)%m;
	if (j<m*m){
            temp=A[j]*X[k];
            atomic_Add(&AA[l], temp); 
	}
    }
}

__global__ void kernel_zero(double *AA, int m) {
    int i,j;
    int t = blockDim.x*gridDim.x;
    int id = threadIdx.x + blockIdx.x*blockDim.x;

    for (i=0;i<m/t;i++){
        j=id+i*t;
        AA[j]=0.0;
    }
}

void cgsolver( double *A, double *b, double *x, int m, int n ){
        double * d_r;
        double * d_rt;
        double * d_p;
        double * d_pt;
        double * d_Ap;
        double * d_tmp;

        double alpha, temp;
	double rsold, rsnew;
        int incx = 1;
        int incy = 1;
        int lda = m;
        double al = 1.;
        double be = 0.;

        int k = 0;

        double *d_A, *d_x, *d_b;

	cudaMalloc((void**)&d_r  ,sizeof(double)*n);
        cudaMalloc((void**)&d_rt ,sizeof(double)*n);
        cudaMalloc((void**)&d_p  ,sizeof(double)*n);
        cudaMalloc((void**)&d_pt ,sizeof(double)*n);
        cudaMalloc((void**)&d_Ap ,sizeof(double)*n);
        cudaMalloc((void**)&d_tmp,sizeof(double)*n);
	cudaMalloc((void**)&d_b  ,sizeof(double)*n);
        cudaMalloc((void**)&d_A  ,sizeof(double)*n*m);
        cudaMalloc((void**)&d_x  ,sizeof(double)*n);

        cudaMemcpy(d_A, A, sizeof(double)*n*m, cudaMemcpyHostToDevice);
        cudaMemcpy(d_x, x, sizeof(double)*n,   cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(double)*n,   cudaMemcpyHostToDevice);

        cublasHandle_t cublasHandle = 0;
        cublasCreate_v2(&cublasHandle);

	kernel_mv<<<10000,1000>>>(d_A, d_x, d_Ap, m); cudaDeviceSynchronize();
	temp=-1.0;

	cublasDcopy_v2(cublasHandle, n, d_b, 1, d_tmp, 1); cudaDeviceSynchronize();
        cublasDaxpy_v2(cublasHandle, n, &temp, d_Ap, 1, d_tmp, 1); cudaDeviceSynchronize();
        cublasDcopy_v2(cublasHandle, n, d_tmp, 1, d_r, 1); cudaDeviceSynchronize();
        cublasDcopy_v2(cublasHandle, n, d_r, 1, d_p, 1); cudaDeviceSynchronize();
        cublasDcopy_v2(cublasHandle, n, d_r, 1, d_rt, 1); cudaDeviceSynchronize();
        cublasDcopy_v2(cublasHandle, n, d_p, 1, d_pt, 1); cudaDeviceSynchronize();

	cublasDdot_v2(cublasHandle, n, d_r, incx, d_rt, incy, &rsold); cudaDeviceSynchronize();

        while ( k < n ){
		kernel_zero<<<10,1000>>>(d_Ap, m); cudaDeviceSynchronize();
		kernel_mv<<<10000,1000>>>(d_A, d_p, d_Ap, m); cudaDeviceSynchronize();
		
		cublasDdot_v2(cublasHandle, n, d_pt, incx, d_Ap, incy, &temp); cudaDeviceSynchronize();
		alpha=rsold/fmax(temp, NEARZERO);

                cublasDaxpy_v2(cublasHandle, n, &alpha, d_p, 1, d_x, 1); cudaDeviceSynchronize();
		temp = -alpha;
                cublasDaxpy_v2(cublasHandle, n, &temp, d_Ap, 1, d_r, 1); cudaDeviceSynchronize();
                cublasDdot_v2(cublasHandle, n, d_r, incx, d_r, incy, &rsnew); cudaDeviceSynchronize();

                if ( sqrt(rsnew) < TOLERANCE ) break;

                cublasDcopy_v2(cublasHandle, n, d_r, 1, d_tmp, 1); cudaDeviceSynchronize();
		temp = rsnew/rsold;
                cublasDaxpy_v2(cublasHandle, n, &temp, d_p, 1, d_tmp, 1); cudaDeviceSynchronize();
                cublasDcopy_v2(cublasHandle, n, d_tmp, 1, d_p, 1); cudaDeviceSynchronize();
                cublasDcopy_v2(cublasHandle, n, d_p, 1, d_pt, 1); cudaDeviceSynchronize();
		rsold = rsnew;

                k++;
        }

        printf("\t[STEP %d] residual = %E\n",k,sqrt(rsold));

	cudaFree(d_A);
	cudaFree(d_x);
        cudaFree(d_r);
        cudaFree(d_rt);
        cudaFree(d_p);
        cudaFree(d_pt);
        cudaFree(d_Ap);
        cudaFree(d_tmp);
	cublasDestroy(cublasHandle);
}

/*
Sparse version of the cg solver
*/

extern "C" void cgsolver_sparse( double *Aval, int *Irn, int *Jcn, double *b, double *x, int n){
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

extern "C" void smvm(int m, const double* val, const int* col, const int* row, const double* x, double* y)
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

extern "C" double * init_source_term(int n, double h){
	double * f;
	int i;
	f  = (double*) malloc(n*sizeof(double*));

	for(i = 0; i < n; i++) {
		f[i] = (double)i * -2. * M_PI * M_PI * sin(10.*M_PI*i*h) * sin(10.*M_PI*i*h);
	}
	return f;
}
