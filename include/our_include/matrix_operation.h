/*
 * matrix_operation.h
 *
 *  Created on: Mar 27, 2019
 *      Author: yll
 */

#ifndef MATRIX_OPERATION_H_
#define MATRIX_OPERATION_H_

void matrix_LU_decomp(double *l, double *u, double *a, int n, int prin);

void matrix_inverse_LU(double *a, int n);

#endif /* MATRIX_OPERATION_H_ */
