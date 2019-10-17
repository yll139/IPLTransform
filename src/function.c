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
	for (int i = 0; i < size; ++i) {
        for (int j = i+1; j < size; ++j) {
        	X[i] += A[size*i+j] * X[j];
        }
    }
	for (int i = 0; i < size; ++i){
		X[i] *= A[size*i+i];
	}

	for (int i = size-1; i >= 0; --i){
		for (int j = i-1; j >= 0; --j)
			X[i] += A[size*i+j] * X[j];
	}
}

void OUR_init_opt(double *A, int size) {
	double *D = (double *)malloc(sizeof(double) * size*size);

	for (int i = 0; i < size*size; ++i)
		D[i] = 0;

	matrix_inverse_LU(A, size);

	for (int i = 0; i < size; ++i) {
		double tmp = A[i*size+i];

		for (int k = 0; k < size; ++k) {
			A[i*size+k] /= tmp;
			D[i*size+k] = 1 / tmp;
		}
		D[i*size+i] = 1 - 1 / tmp;

		for (int j = 0; j < size; ++j) {
			if (j != i) {
				double tmp = A[j*size+i];
				for (int k = 0; k < size; ++k) {
					A[j*size+k] -= A[i*size+k] * tmp;
				}
				D[i*size+j] *= tmp;
			}
		}
	}

	for (int i = 0; i < size*size; ++i)
		A[i] = D[i];

	free(D);
}

void OUR_init(double *A, int size) {
	double *para = (double *)malloc(sizeof(double) * size*size);

	for (int i=size-1, s=0; i > 0; --i, ++s) {
		double *B = (double *)malloc(sizeof(double) * i*i);
		double *C = (double *)malloc(sizeof(double) * i*i);

		//基矩阵B, 副本C
		for (int j = 0; j < i; ++j) {
			for (int k = 0; k < i; ++k) {
				int row = s + 1 + j;
				int col = s + 1 + k;
				B[j*i+k] = A[row*size+col];
				C[j*i+k] = A[row*size+col];
			}
		}

		//分母(基矩阵的行列式)=down
		double *l = (double *)malloc(sizeof(double) * i*i);
		double *u = (double *)malloc(sizeof(double) * i*i);
		double down = 1.0;
		matrix_LU_decomp(l, u, B, i, 0);
		for (int j = 0; j < i; ++j) {
			down *= u[j*i+j];
		}

		//需组合得到的列向量D
		double *D = (double *)malloc(sizeof(double) * i);
		int col = s;
		for (int j = 0; j < i; ++j) {
			int row = s + 1 + j;
			D[j] = A[row*size+col];
		}

		//得到参数矩阵para的右上半部分
		for (int j = 0; j < i; ++j) {
			//替换B矩阵的第j列,为列向量D
			for (int k = 0; k < i; ++k) {
				B[k*i+j] = D[k];
			}

			double upper = 1.0;
			matrix_LU_decomp(l, u, B, i, 0);
			for (int j = 0; j < i; ++j) {
				upper *= u[j*i+j];
			}

			int tmp = s*size + s + 1 + j;
			para[tmp] = upper / down;

			//由副本C恢复基矩阵B
			for (int k = 0; k < i*i; ++k) {
				B[k] = C[k];
			}
		}

		//得到参数矩阵para的左下半部分
		//para[s*size+s] = A[s*size+s];
		para[s*size+s] = A[s*size+s] - 1;
		for (int j = 0; j < i; ++j) {
			int tmp = s*size + s + 1 + j;
			para[s*size+s] -= A[tmp] * para[tmp];
		}

		for (int k = s-1; k >= 0; --k) {
			double tmp = 0.0;
			for (int j = 0; j < i; ++j) {
				int tmp1 = k*size + s + 1 + j;
				int tmp2 = s*size + s + 1 + j;
				tmp += A[tmp1] * para[tmp2];
			}
			para[s*size+k] = (A[k*size + s] - tmp);
		}

		//释放内存
		free(B);
		free(C);
		free(D);
		free(l);
		free(u);
	}

	//参数的最后一列（由全单元素向量构成)
	int tmp = (size - 1) * size;
	para[tmp+size-1] = A[tmp+size-1] - 1;

	for (int i = size-2; i>=0; --i) {
		para[tmp+i] = A[i*size+size-1];
	}

	//赋值回矩阵A
	for (int i = 0; i < size*size; ++i)
		A[i] = -para[i];

	//释放内存
	free(para);
}

void OUR_function(double *X, double *A, int size) {
	for (int i = 0; i < size; ++i) {
		double tmp = X[i];
		for (int j=0; j < size; ++j) {
			X[j] -= A[i*size+j] * tmp;
		}
	}
}
