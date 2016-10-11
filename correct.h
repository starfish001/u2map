#ifndef _CORRECT_H_
#define _CORRECT_H_

int readline2(FILE * fp,char *mark,char *mark_name,int SAMP,int grouptype);
double mle2kosambi(double r);
double kosambi2mle(double k);
int locate(double distAX,double distBX,double distAB);
double new_MLE(double distAX,double distAY,double distBX,double distBY,double distAB);
void correct_CP_MLE(char *markFile,char *inputPrefix,double accuracy);

#endif
