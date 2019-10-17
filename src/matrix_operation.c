/*
 * matrix_operation.c
 *
 *  Created on: Mar 25, 2019
 *      Author: yll
 */
#include <stdio.h>
#include <stdlib.h>

void matrix_LU_decomp(double *l, double *u, double *a, int n, int prin) {
    int i, j, k;
    double s;

    for (i = 0; i < n*n; i++)
    	l[i] = u[i] = 0;

    for (i = 0; i < n;i++)      //计算l矩阵对角线
        l[i*n+i] = 1;

    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            s = 0;
            for (k = 0; k < i; k++) {
                s += l[i*n+k] * u[k*n+j];
            }
            u[i*n+j] = a[i*n+j] - s;      //按行计算u值
        }

        for (j = i + 1; j < n; j++) {
            s = 0;
            for (k = 0; k < i; k++) {
                s += l[j*n+k] * u[k*n+i];
            }
            l[j*n+i] = (a[j*n+i] - s) / u[i*n+i];      //按列计算l值
        }
    }

    if (prin) {
    	printf("matrix L:\n");
    	for (i = 0; i < n; i++) {
    		for (j = 0; j < n; j++) {
    			printf("%.2f\t", l[i*n+j]);
    		}
    		printf("\n");
    	}

    	printf("matrix U:\n");
    	for (i = 0; i < n; i++) {
    		for (j = 0; j < n; j++) {
    			printf("%.2f\t", u[i*n+j]);
    		}
    		printf("\n");
    	}
    }
}

void matrix_inverse_LU(double *a, int n) {
    double *l, *u;
    double *l_inverse, *u_inverse;
    double *a_inverse;
    int i, j, k;
    double s;

    l = (double *)malloc(sizeof(double) * n*n);
    u = (double *)malloc(sizeof(double) * n*n);
    l_inverse = (double *)malloc(sizeof(double) * n*n);
    u_inverse = (double *)malloc(sizeof(double) * n*n);
    a_inverse = (double *)malloc(sizeof(double) * n*n);

    for (i = 0; i < n*n; i++) {
    	l_inverse[i] = u_inverse[i] = 0;
    	a_inverse[i] = 0;
    }

    matrix_LU_decomp(l, u, a, n, 0);


    for (i = 0; i < n; i++) {        //按行序，行内从高到低，计算l的逆矩阵
        l_inverse[i*n+i] = 1;
    }
    for (i= 1; i < n; i++) {
        for (j = 0; j < i; j++) {
            s = 0;
            for (k = 0; k < i; k++) {
                s += l[i*n+k] * l_inverse[k*n+j];
            }
            l_inverse[i*n+j] = -s;
        }
    }

    for (i = 0; i < n; i++) {                  //按列序，列内按照从下到上，计算u的逆矩阵
        u_inverse[i*n+i] = 1 / u[i*n+i];
    }
    for (i = 1; i < n; i++) {
        for (j = i - 1; j >= 0; j--) {
            s = 0;
            for (k = j + 1 ; k <= i; k++) {
                s += u[j*n+k] * u_inverse[k*n+i];
            }
            u_inverse[j*n+i] = -s / u[j*n+j];
        }
    }

    for (i = 0;i < n;i++) {           //计算矩阵a的逆矩阵
        for (j = 0;j < n;j++) {
            for (k = 0;k < n;k++) {
                a_inverse[i*n+j] += u_inverse[i*n+k] * l_inverse[k*n+j];
            }
        }
    }

    for (i = 0; i < n; i++)
    	for (j = 0; j < n; j++)
    		a[i*n+j] = a_inverse[i*n+j];

    free(l);
    free(u);
    free(l_inverse);
    free(u_inverse);
    free(a_inverse);
}
