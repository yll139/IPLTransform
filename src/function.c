/*
 * function.c
 *
 *  Created on: Mar 25, 2019
 *      Author: yll
 */
#include <stdio.h>
#include <stdlib.h>
#include "matrix_operation.h"

void LDUBY_init(double *A, int size) {
	double *l = (double *)malloc(sizeof(double) * size*size);
	double *u = (double *)malloc(sizeof(double) * size*size);

	matrix_LU_decomp(l, u, A, size, 0);

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			A[i*size+j] = (i>j) ? l[i*size+j] : u[i*size+j]/u[i*size+i];
			if (j == i)
				A[i*size+i] = u[i*size+i];
		}
	}
	free(l);
	free(u);
}

void LUBY_function(double *X, double *A, int size) {
	double memory_one;

	for (int i = 0; i < size; ++i) {
        for (int j = i+1; j < size; ++j) {
        	memory_one = A[size*i+j] * X[j];
        	X[i] += memory_one;
        }
    }

	for (int i = 0; i < size; ++i)
		X[i] *= A[size*i+i];

	for (int i = size-1; i >= 0; --i) {
		for (int j = i-1; j >= 0; --j) {
			memory_one = A[size*i+j] * X[j];
			X[i] += memory_one;
		}
	}

}

void OUR_init_opt(double *A, int size) {
	double *D = (double *)malloc(sizeof(double) * size*size);

	for (int i = 0; i < size*size; ++i)
		D[i] = 0;

	matrix_inverse_LU(A, size);

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			if (i == j)	continue;

			if (A[j*size+i]) {
				D[i*size+j] = -A[i*size+i] / A[j*size+i];
				for (int k = 0; k < size; ++k)
					A[j*size+k] = A[j*size+k] * D[i*size+j] + A[i*size+k];
			}
		}

		D[i*size+i] = 1 / A[i*size+i];
		for (int j = 0; j < i; ++j)
			D[j*size+j] /= D[i*size+j];
	}

	for (int i = 0; i < size*size; ++i)
		A[i] = D[i];

	free(D);
}

void OUR_function(double *X, double *A, int size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < i; ++j) {
			X[j] *= A[i*size+j];
			X[j] += X[i];
		}

		for (int j = i+1; j < size; ++j) {
			X[j] *= A[i*size+j];
			X[j] += X[i];
		}
	}

	for (int i = 0; i < size; ++i)
		X[i] *= A[i*size+i];
}
