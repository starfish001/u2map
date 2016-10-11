#ifndef _PREPARE_H_
#define _PREPARE_H_

int step1_prepare(char *inputFileName,char *outputFileName,int grouptype,float accuracy);
void calculate(double accuracy,double r,float eachN[],double *result,int grouptype,int marktype1,int marktype2);
int readline(FILE * fp,char *mark,char *mark_name,int SAMP,float *percentage,int grouptype);
double calLOD(double lrt,double mle);
void statMark(char *make_name1,char *make_name2,char *mark1,char *mark2,float *eachN,int SAMP,float percentage1[],float percentage2[],int grouptype);
void errorPrint();

#endif
