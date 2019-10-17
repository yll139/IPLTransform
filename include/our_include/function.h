/*
 * function.h
 *
 *  Created on: Mar 27, 2019
 *      Author: yll
 */

#ifndef FUNCTION_H_
#define FUNCTION_H_

void LDUBY_init(double *A, int size);

void LUBY_function(double *X, double *A, int size);

void OUR_init(double *A, int size);

void OUR_init_opt(double *A, int size);

void OUR_function(double *X, double *A, int size);

#endif /* FUNCTION_H_ */
