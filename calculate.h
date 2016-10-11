#ifndef _CALCULATE_H_

#define _CALCULATE_H_
double corrcteMLE(double r);

double MLE1(float n2,float n3,double n);
double VAR1(float n1,float n2,float n3,float n4,double n);
double LRT1(float n1,float n2,float n3,float n4,double n,double r);
void f1dh(float n[],double *result);

double MLE2(double n1,double n2,double n3,double n4,double n);
double VAR2(double n2,double n3,double n,double r);
double LRT2(double n2,double n3,double n,double r);
void f1ril(float n[],double *result);

double MLE3(double n1,double n2,double n3,double n4,double n,double r0,double accuracy);
double VAR3(double n1,double n2,double n3,double n4,double n,double r);
double LRT3(double n1,double n2,double n3,double n4,double n,double r);
void p1bc1ril(float n[],double *result,double r,double accuracy);

double MLE4(double n1,double n2,double n3,double n4,double n,double r0,double accuracy);
double VAR4(double n1,double n2,double n3,double n4,double n,double r);
double LRT4(double n1,double n2,double n3,double n4,double n,double r);
void p2bc1ril(float n[],double *result,double r,double accuracy);

double MLE5(double n1,double n2,double n3,double n4,double n,double r0,double accuracy);
double VAR5(double n1,double n2,double n3,double n4,double n,double r);
double LRT5(double n1,double n2,double n3,double n4,double n,double R);
void p1bc2ril(float n[],double *result,double r,double accuracy);

double MLE6(double n1,double n2,double n3,double n4,double n,double r0,double accuracy);
double VAR6(double n1,double n2,double n3,double n4,double n,double r);
double LRT6(double n1,double n2,double n3,double n4,double n,double R);
void p2bc2ril(float n[],double *result,double r,double accuracy);

double MLE7(double n1,double n2,double n3,double n4,double n,double r0,double accuracy);
double VAR7(double n1,double n2,double n3,double n4,double n,double r);
double LRT7(double n1,double n2,double n3,double n4,double n,double R);
void p1bc1dh(float n[],double *result,double r,double accuracy);

double MLE8(double n1,double n2,double n3,double n4,double n,double r0,double accuracy);
double VAR8(double n1,double n2,double n3,double n4,double n,double r);
double LRT8(double n1,double n2,double n3,double n4,double n,double R);
void p2bc1dh(float n[],double *result,double r,double accuracy);

double MLE9(double n1,double n2,double n3,double n4,double n,double r0,double accuracy);
double VAR9(double n1,double n2,double n3,double n4,double n,double r);
double LRT9(double n1,double n2,double n3,double n4,double n,double R);
void p1bc2dh(float n[],double *result,double r,double accuracy);

double MLE10(double n1,double n2,double n3,double n4,double n,double r0,double accuracy);
double VAR10(double n1,double n2,double n3,double n4,double n,double r);
double LRT10(double n1,double n2,double n3,double n4,double n,double R);
void p2bc2dh(float n[],double *result,double r,double accuracy);

double MLE11(double n2,double n4,double n);
double VAR11(double n1,double n2,double n4,double n5,double r);
double LRT11(double n1,double n2,double n4,double n5,double n,double r);
void p1bc1f1(float n[],double *result,double r,double accuracy);

double MLE12(double n6,double n8,double n);
double VAR12(double n5,double n6,double n8,double n9,double n,double r);
double LRT12(double n5,double n6,double n8,double n9,double n,double r);
void p2bc1f1(float n[],double *result,double r,double accuracy);

double MLE13(double n5,double N1,double N2,double r,double accuracy);
double VAR13(double n5,double N1,double N2,double r);
double LRT13(double n5,double N1,double N2,double n,double r);
void f2(float n[],double *result,double r,double accuracy);

void f3g(float n[],double *result,double R,double accuracy);

double MLE15(double n1,double n5,double n);
double VAR15(double n1,double n2,double n4,double n5,double n,double r);
double LRT15(double n1,double n2,double n4,double n5,double R);
void p1bc2f1(float n[],double *result,double r,double accuracy);

double MLE16(double n5,double n9,double n);
double VAR16(double n5,double n6,double n8,double n9,double n,double r);
double LRT16(double n5,double n6,double n8,double n9,double R);
void p2bc2f1(float n[],double *result,double r,double accuracy);

void p1bc1f2(float n[],double *result,double R,double accuracy);

void p2bc1f2(float n[],double *result,double R,double accuracy);

void p1bc2f2(float n[],double *result,double R,double accuracy);

void p2bc2f2(float n[],double *result,double R,double accuracy);

double callLOD(double lrt,double mle);

void get_mi_matrix(int phenotype,char mark,double (*mi)[1]);
void get_H(double r,double *H_point,int phase);
void get_D(double *D_point,int phase);
int get_I_row(int phenotype);
void get_I_array(double *I_point,int row,int phenotype);
void cp(int phenotype1,int phenotype2,char *mark1_array,char *mark2_array,double *result,double r0,double accuracy,int SAMP,int LOD_flag);

#endif
