/*
 * perf.c
 *
 *  Created on: Mar 26, 2019
 *      Author: yll
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include "mkl.h"

#include "function.h"

//#define DEBUG
#define MKL_FUN

int main(int argc,  char *argv[])
{
    double *Y;
    double *M, *M_bak;
    double *X, *X_bak;
    double s_initial, s_elapsed0, s_elapsed2, s_elapsed3;

    double max_error, ave_error;
	
	int Size = (argc > 1) ? atoi(argv[1]) : 50;//the default size of square matrix M is 50x50.
    int LOOP_COUNT = (argc > 2) ? atoi(argv[2]) : 10000;//the default size of LOOP_COUNT is 10000.
	
    printf("Y = M x X, \nwhere the size of square matrix M is %dx%d and the size of vector X is %dx1.\n\n", 
           Size, Size, Size);

#ifdef MKL_FUN
    double s_elapsed1;
#endif
    
#ifdef DEBUG
    Size = 2;
    LOOP_COUNT = 1;
#endif

//++++++++++++++++++++++++++++++++++++Memory++++++++++++++++++++++++++++++++++++
    M       =  (double *)mkl_malloc(sizeof(double) * Size*Size, 64);        //++
    M_bak   =  (double *)mkl_malloc(sizeof(double) * Size*Size, 64);        //++
    X       =  (double *)mkl_malloc(sizeof(double) * Size, 64);             //++
    X_bak   =  (double *)mkl_malloc(sizeof(double) * Size, 64);             //++
    Y       =  (double *)mkl_malloc(sizeof(double) * Size, 64);             //++
    if (M == NULL || X == NULL || Y == NULL) {                              //++
        printf( "\nERROR\n\n");                                             //++
        mkl_free(M);                                                        //++
        mkl_free(M_bak);                                                    //++
        mkl_free(X);                                                        //++
        mkl_free(X_bak);                                                    //++
        mkl_free(Y);                                                        //++
        return 1;                                                           //++
    }                                                                       //++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//++


//+++++++++++++++++++++++++++++Intializing matrix data++++++++++++++++++++++//++
	srand((unsigned)time(NULL));                              				//++
    for (int i = 0; i < (Size*Size); i++) {                                 //++
        M[i] = ((double)rand()) / RAND_MAX;                                 //++
        M_bak[i] = M[i];                                                    //++
    }                                                                       //++
                                                                            //++
    for (int i = 0; i < Size; i++) {                                        //++
        X[i] = ((double)rand()) / RAND_MAX;                                 //++
        X_bak[i] = X[i];                                                    //++
        Y[i] = 0.0;                                                         //++
    }                                                                       //++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//++

#ifdef DEBUG
    printf("==========M*X==========\n");
    for (int i = 0; i < Size; i++){
        printf("[");
        for (int j = 0; j < Size; j++) {
            printf("%.2f\t", M[i*Size+j]);
        }
        printf("] [%.2f]\n", X[i]);
    }
    printf("=======================\n\n");
#endif

//+++++++++++++++++++++++++++++++triple nested loop+++++++++++++++++++++++++++++++++++++++++//++
    printf ("1. Measuring performance of matrix product using triple nested loop \n");      //++
    s_elapsed0 = 0;                                                                         //++
    for (int l = 0; l < LOOP_COUNT; ++l) {                                                  //++
        for (int i = 0; i < Size; i++)                                                      //++
            Y[i] = 0.0;                                                                     //++
                                                                                            //++
        s_initial = dsecnd();                                                               //++
        for (int i = 0; i < Size; i++) {                                                 	//++
            for (int j = 0; j < Size; j++) {                                             	//++
                Y[i] += M[i*Size+j] * X[j];                                                 //++
            }                                                                               //++
        }                                                                                   //++
        s_elapsed0 += dsecnd() - s_initial;                                                 //++
    }                                                                                       //++
    printf ("** Matrix multiplication using triple nested loop completed**\n"               //++
            "=> at %.5f microseconds <=\n\n", 1000000 * s_elapsed0/LOOP_COUNT);             //++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//++


#ifdef DEBUG
    printf("===Y===\n");
    for (int i = 0; i < Size; i++){
        printf("[%.2f]\n", Y[i]);
    }
    printf("=======\n\n");
#endif


#ifdef MKL_FUN
//+++++++++++++++++++++++++++++++Intel(R) MKL dgemm function++++++++++++++++++++++++++++++++++//++
    printf ("2. Measuring performance of matrix product using Intel(R) MKL dgemm function \n" //++
            " via CBLAS interface \n\n"); 													  //++
	mkl_set_num_threads(1); 																  //++
	s_elapsed1 = 0; 																		  //++
    for (int i = 0; i < LOOP_COUNT; i++) { 													  //++
        for (int i = 0; i < Size; i++)                                                        //++
            Y[i] = 0.0;                                                       				  //++
        																					  //++
        s_initial = dsecnd();                                                       		  //++
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,                                //++
                    Size, 1, Size, 1.0, M, Size, X, 1, 0.0, Y, 1);                            //++
        s_elapsed1 += dsecnd() - s_initial;                                                   //++
        																					  //++
    }                                                       								  //++
    printf ("** Matrix multiplication using Intel(R) MKL dgemm completed**\n"                 //++
            "=> at %.5f microseconds <=\n\n", 1000000 * s_elapsed1/LOOP_COUNT);               //++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//++
#endif

    printf("-------------Initializing-----------\n");
    LDUBY_init(M, Size);
    printf("---------Initialization completed-----\n");

#ifdef DEBUG
    printf("==========LDU*X==========\n");
    for (int i = 0; i < Size; i++){
        printf("[");
        for (int j = 0; j < Size; j++) {
            printf("%.2f\t", M[i*Size+j]);
        }
        printf("] [%.2f]\n", X[i]);
    }
    printf("=======================\n");
#endif

//++++++++++++++++++++++++++++++++++++LUBY's method+++++++++++++++++++++++++++++++++++++++++//++
    printf ("3. Measuring performance of matrix-vector product using LUBY's method \n");	//++
    s_elapsed2 = 0;                                                                         //++
    for (int i = 0; i < LOOP_COUNT; i++) {                                                  //++
                                                                                            //++
        for (int j = 0; j < Size; j++)                                                      //++
            X[j] = X_bak[j];                                                                //++
                                                                                            //++
        s_initial = dsecnd();                                                               //++
        LUBY_function(X, M, Size);                                                          //++
        s_elapsed2 += dsecnd() - s_initial;                                                 //++
    }                                                                                       //++
    printf ("** Matrix multiplication using LUBY's method completed**\n"                    //++
            " => at %.5f microseconds <=\n\n", 1000000 * s_elapsed2/LOOP_COUNT);            //++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//++


#ifdef DEBUG
    printf("===Y===\n");
    for (int i = 0; i < Size; i++){
        printf("[%.2f]\n", X[i]);
    }
    printf("=======\n\n");
#endif

    //Check
    for (int i = 0; i < Size; i++) {
        double tmp = fabs(X[i] - Y[i]);
        if (tmp > 0.00001){
            printf("Fail\n");
            return 2;
        }
    }
	max_error = 0.0;
    ave_error = 0.0;
    for (int i = 0; i < Size; i++) {
        double tmp = fabs( (X[i] - Y[i]) / Y[i] );
        if (tmp > max_error){
        	max_error = tmp;
        }
        ave_error += tmp;
    }
    printf("Max error: %g Average error: %g\n", max_error, ave_error / Size);


    for (int j = 0; j < Size; j++)
        X[j] = X_bak[j];

    for (int j = 0; j < Size*Size; j++)
        M[j] = M_bak[j];

    printf("-------------Initializing-----------\n");
    OUR_init_opt(M, Size);
    printf("---------Initialization completed-----\n");

#ifdef DEBUG
    printf("==========OUR*X==========\n");
    for (int i = 0; i < Size; i++){
        printf("[");
        for (int j = 0; j < Size; j++) {
            printf("%.2f\t", M[i*Size+j]);
        }
        printf("] [%.2f]\n", X[i]);
    }
    printf("=======================\n\n");
#endif

//++++++++++++++++++++++++++++++++++++Our algorithm+++++++++++++++++++++++++++++++++++++++++//++
    printf ("4. Measuring performance of matrix product using our algorithm \n");			//++
    s_elapsed3 = 0;                                                                         //++
    for (int i = 0; i < LOOP_COUNT; i++) {                                                  //++
                                                                                            //++
        for (int j = 0; j < Size; j++)                                                      //++
            X[j] = X_bak[j];                                                                //++
                                                                                            //++
        s_initial = dsecnd();                                                               //++
        OUR_function(X, M, Size);                                                           //++
        s_elapsed3 += dsecnd() - s_initial;                                                 //++
    }                                                                                       //++
    printf ("** Matrix multiplication using our algorithm completed**\n"                    //++
            "=> at %.5f microseconds <=\n\n", 1000000 * s_elapsed3/LOOP_COUNT);             //++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//++


#ifdef DEBUG
    printf("===Y===\n");
    for (int i = 0; i < Size; i++){
        printf("[%.2f]\n", X[i]);
    }
    printf("=======\n\n");
#endif

    //Check
    for (int i = 0; i < Size; i++) {
        double tmp = fabs(X[i] - Y[i]);
        if (tmp > 0.00001){
            printf("Fail\n");
            return 3;
        }
    }
    max_error = 0.0;
    ave_error = 0.0;
    for (int i = 0; i < Size; i++) {
        double tmp = fabs( (X[i] - Y[i]) / Y[i] );
        if (tmp > max_error){
        	max_error = tmp;
        }
        ave_error += tmp;
    }
    printf("Max error: %g Average error: %g\n", max_error, ave_error / Size);

    printf("\n===========RATE: %.5f %%==============\n", 100 *((double)(s_elapsed2 - s_elapsed3))/s_elapsed2);



    mkl_free(M);
    mkl_free(M_bak);
    mkl_free(X);
    mkl_free(X_bak);
    mkl_free(Y);

    printf ("Example completed. \n\n");
    return 0;
}
