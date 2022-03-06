#include <iostream>
#include "mkl_lapacke.h"


/* Auxiliary routines prototypes */
extern void print_matrix(char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda);


int main(int argc, char * argv[]){
	// 1. Basic linear solving, Ax = b
	/*
	A = 1 2 3
	    4 5 6
		7 8 10

	b = 3
	    3 
		4

	x = -2
	    1
		1
	*/
	MKL_INT n = 3, nrhs = 1, lda = 3, ldb = 3, info;
	MKL_INT ipiv[3];
	double a[9] = {
		1, 4, 7, 2, 5, 8, 3, 6, 10
	};
	double b[3] = { 3, 3, 4 };
	printf("LAPACKE_dgesv (column-major, high-level) Example Program Results\n");
	/* Solve the equations A*X = B */
	// dgesv ref: https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/lapack-routines/lapack-linear-equation-routines/lapack-linear-equation-driver-routines/gesv.html 
	info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
	
	/* Check for the positive definiteness */
	if (info > 0) {
		printf("The leading minor of order %i is not positive ", info);
		printf("definite;\nthe solution could not be computed.\n");
		return -1;
	}
	/* Print solution */
	print_matrix("Solution", n, nrhs, b, ldb);
	/* Print details of Cholesky factorization */
	print_matrix("Details of Cholesky factorization", n, n, a, lda);

	return 0;
}


/* Auxiliary routine: printing a matrix */
void print_matrix(char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda) {
	MKL_INT i, j;
	printf("\n %s\n", desc);
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) printf(" %6.2f", a[i + j * lda]);
		printf("\n");
	}
}
