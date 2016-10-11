#ifndef _MATRIX_H_
#define _MATRIX_H_

void matrix_elementwise(double *M1,int M1_h,int M1_l,double *M2,int M2_h,int M2_l,double *M3,int M3_h,int M3_l);//an elementwise product of two matrices
void matrix_multi(double *M1,int M1_h,int M1_l,double *M2,int M2_h,int M2_l,double *M3,int M3_h,int M3_l);//矩阵相乘
void matrix_transp(double *M1,int M1_h,int M1_l,double *M2,int M2_h,int M2_l);//矩阵转置
void matrix_printf(double *M1,int M1_h,int M1_l);//矩阵输出

#endif
