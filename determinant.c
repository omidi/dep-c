#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 

double determinant(double *A, int n) {
    // will be holding the Cholesky decomposition of the matrix A
    double *L = (double*)calloc(n * n, sizeof(double)); 
    if (L == NULL)
        exit(EXIT_FAILURE);

    double det = 0.; // will hold the determinant of the matrix 
    int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j < (i+1); j++) {
            double s = 0;
            for (k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
                           sqrt(A[i * n + i] - s) :
                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
	    if(i==j){
	      det += 2.0 * log(L[i * n + j]);
	    }
        }
    free(L);
    return det;
}


void cfun(const void * indatav, size_t size, void * outdatav) {
    //void cfun(const double * indata, int rowcount, int colcount, double * outdata) {
    size_t i;
    const double * indata = (double *) indatav;
    double * outdata = (double *) outdatav;
    puts("Here we go!");
    for (i = 0; i < size; ++i) {
        outdata[i] = indata[i] * 2.0;
    }
    puts("Done!");
}
