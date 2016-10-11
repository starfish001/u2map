#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

//an elementwise product of two matrices
void matrix_elementwise(double *M1,int M1_h,int M1_l,double *M2,int M2_h,int M2_l,double *M3,int M3_h,int M3_l){
        if(M1_h != M2_h || M1_l != M2_l){fprintf(stderr,"Error：matrix 1 row not equal matrix 2 line.\n");exit(EXIT_FAILURE);}
        if(M3_h != M1_h || M3_l != M1_l){fprintf(stderr,"Error: matrix 3 row or line error.\n");exit(EXIT_FAILURE);}

        double (*A)[M1_l] = (double (*)[M1_l])M1;
        double (*B)[M2_l] = (double (*)[M2_l])M2;
        double (*C)[M3_l] = (double (*)[M3_l])M3; //C= A o B

        int i,j;
        for(i=0;i<M1_h;i++){
                for(j=0;j<M2_l;j++){
                       (*(*(C+i)+j)) = (*(*(A+i)+j)) * (*(*(B+i)+j));
                }
        }
}


//矩阵相乘
void matrix_multi(double *M1,int M1_h,int M1_l,double *M2,int M2_h,int M2_l,double *M3,int M3_h,int M3_l){
	if(M1_l != M2_h){fprintf(stderr,"Error：matrix 1 row not equal matrix 2 line.\n");exit(EXIT_FAILURE);}
	if(M3_h != M1_h || M3_l != M2_l){fprintf(stderr,"Error: matrix 3 row or line error.\n");exit(EXIT_FAILURE);}

	double (*A)[M1_l] = (double (*)[M1_l])M1;
	double (*B)[M2_l] = (double (*)[M2_l])M2;
	double (*C)[M3_l] = (double (*)[M3_l])M3; //C= A*B C的列数等于B的列数

	int i,j,x;
	for(i=0;i<M1_h;i++){
		for(j=0;j<M2_l;j++){
			double temp = 0;
			for(x=0;x<M1_l;x++){
				temp += (*(*(A+i)+x)) * (*(*(B+x)+j));
			}
			*(*(C+i)+j) = temp;
		}
	}
}

//矩阵转置
void matrix_transp(double *M1,int M1_h,int M1_l,double *M2,int M2_h,int M2_l){
	if(M1_h != M2_l || M1_l != M2_h){fprintf(stderr,"Error: matrix 2 row or line error.\n");exit(EXIT_FAILURE);}
	double (*A)[M1_l] = (double (*)[M1_l])M1;
	double (*B)[M2_l] = (double (*)[M2_l])M2;
	int i,j;
	for(i=0;i<M2_h;i++){
		for(j=0;j<M2_l;j++){
			*(*(B+i)+j) = *(*(A+j)+i);
		}
	}
}

//矩阵输出
void matrix_printf(double *M1,int M1_h,int M1_l){
	double (*A)[M1_l] = (double (*)[M1_l])M1;
	int i,j;
	for(i=0;i<M1_h;i++){
		for(j=0;j<M1_l;j++){
			printf("%lf\t",*(*(A+i)+j));
		}
		printf("\n");
	}
}

/*
int main(){
	double A[2][4] = {{1,2,3,4},{1,2,3,4}};
	double B[4][2] = {{1,1},{2,2},{3,3},{4,4}};
	double C[2][2];
	double D[2][4];

	double * M1 = (double *)A;
	double * M2 = (double *)B;
	double * M3 = (double *)C;
	double * M4 = (double *)D;
	matrix_multi(M1,2,4,M2,4,2,M3,2,2);
	matrix_printf(M3,2,2);
	printf("hhh\n");
	matrix_transp(M2,4,2,M4,2,4);
	matrix_printf(M4,2,4);
	return 0;
}*/

