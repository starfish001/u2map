#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "calculate.h"
//20150517 针对计算MLE值出现死循环的问题，界定循环数为1000.取temp值最小的R值。（该判断方式也有问题，给出的MLE未必是最接近的）
//20150521 针对不规范MLE值纠正，（小于0、大于0.5、nan、inf的MLE都纠正为0.5
//This  is for group "F1DH"


//纠正函数
double corrcteMLE(double r)
{
	if(r >= 0 && r <= 0.5){return r;}
	else if(r > 0.5){return 1 - r;}
	else return 0.5;

}

double MLE1(float n2,float n3,double n)
{
	double r = (n2+n3)/n;
	r = corrcteMLE(r);
	return r; 
}

double VAR1(float n1,float n2,float n3,float n4,double n)
{
	return 1/((pow(n,3))/((n1+n4)*(n2+n3)));
}

double LRT1(float n1,float n2,float n3,float n4,double n,double r)
{
	double x,y,z;
	x = pow(0.5,n);
	y = (pow((1-r),(n1+n4)))*(pow(r,(n2+n3)));
	z = log(x/y) * (-2);
	return z;
}

void f1dh(float n[],double *result){
        float n1 = n[0];
        float n2 = n[2];
        float n3 = n[6];
        float n4 = n[8];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE1(n2,n3,all);
        *(result + 1) = VAR1(n1,n2,n3,n4,all);
        *(result + 2) = LRT1(n1,n2,n3,n4,all,*result);
}

//This  is for group "F1RIL"
double MLE2(double n1,double n2,double n3,double n4,double n)
{
        double r = (n2+n3)/(2*(n1+n4));
	r = corrcteMLE(r);
	return r;
}

double VAR2(double n2,double n3,double n,double r)
{
        double x,y,z;
	x = pow(r,2);
	y = pow((1+2*r),2);
	z = (n2+n3)/x - 4*n/y;
	return 1/z;
}

double LRT2(double n2,double n3,double n,double r)
{
        double x,y,z;
        x = pow((1+2*r),n);
	y = pow(2*r,(n2+n3));
	z = log(x/y)*(-2);
	return z;
}

void f1ril(float n[],double *result){
        double n1 = n[0];
        double n2 = n[2];
        double n3 = n[6];
        double n4 = n[8];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE2(n1,n2,n3,n4,all);
        *(result + 1) = VAR2(n2,n3,all,*result);
        *(result + 2) = LRT2(n2,n3,all,*result);
}

//This is for group "P1BC1RIL"
double MLE3(double n1,double n2,double n3,double n4,double n,double r0,double accuracy){
        double d1;//[ln(L(r))]'
	double d2;//[ln(L(r))]''
	double min = 1;
	double minR;
	int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的
	double r = r0;
	double temp;
	do{
		d1 = n1/(1+r) + (n2+n3)/r - n4/(1-r) - (2*n)/(1+2*r);
		d2 = (-1)*n1/(pow((1+r),2)) - (n2*n3)/(pow(r,2)) - n4/(pow((1-r),2)) + 4*n/(pow((1+2*r),2));
		temp = d1/d2;
		r -=temp;
		num++;
		if(fabs(temp) < min){min = fabs(temp);minR = r;}
		if(num > 1000){r = minR;break;}
	}while(fabs(temp)>accuracy);
	r = corrcteMLE(r);
	return r;
}

double VAR3(double n1,double n2,double n3,double n4,double n,double r){
	double z;
	z = n1/(pow((1+r),2)) + (n2*n3)/(pow(r,2)) + n4/(pow((1-r),2)) - 4*n/(pow((1+2*r),2));
	return 1/z;
}

double LRT3(double n1,double n2,double n3,double n4,double n,double r){
        double a,b,c,d,e,z;
	a = n1 * log(3);
	b = n * log(0.5 + r);
	c = n1 * log(1+r);
	d = (n2+n3) * log(r);
	e = n4 * log(1-r);
	z = (-2) * (a+b-c-d-e);
	return z;
}

void p1bc1ril(float n[],double *result,double r,double accuracy){
        double n1 = n[0];
        double n2 = n[2];
        double n3 = n[6];
        double n4 = n[8];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE3(n1,n2,n3,n4,all,r,accuracy);
        *(result + 1) = VAR3(n1,n2,n3,n4,all,*result);
        *(result + 2) = LRT3(n1,n2,n3,n4,all,*result);
}

/*This Headfiel is for group "P2BC1RIL" 
the deference between P2BC1RIL and P1BC1RIL is: n1 <-> n4 
*/
double MLE4(double n1,double n2,double n3,double n4,double n,double r0,double accuracy){
        double d1;//[ln(L(r))]'
	double d2;//[ln(L(r))]''
	double r = r0;
	double temp;
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的

	do{
		d1 = n1/(1+r) + (n2+n3)/r - n4/(1-r) - (2*n)/(1+2*r);
		d2 = (-1)*n1/(pow((1+r),2)) - (n2*n3)/(pow(r,2)) - n4/(pow((1-r),2)) + 4*n/(pow((1+2*r),2));
		temp = d1/d2;
		r -=temp;
//printf("temp=%lf,r=%lf\n",temp,r);
                num++;
                if(fabs(temp) < min){min = fabs(temp);minR = r;}
                if(num > 1000){r = minR;break;}
	}while(fabs(temp)>accuracy);
	r = corrcteMLE(r);
	return r;
}

double VAR4(double n1,double n2,double n3,double n4,double n,double r){
	double z;
	z = n1/(pow((1+r),2)) + (n2*n3)/(pow(r,2)) + n4/(pow((1-r),2)) - 4*n/(pow((1+2*r),2));
	return 1/z;
}

double LRT4(double n1,double n2,double n3,double n4,double n,double r){
        double a,b,c,d,e,z;
	a = n1 * log(3);
	b = n * log(0.5 + r);
	c = n1 * log(1+r);
	d = (n2+n3) * log(r);
	e = n4 * log(1-r);
	z = (-2) * (a+b-c-d-e);
	return z;
}

void p2bc1ril(float n[],double *result,double r,double accuracy){
        double n1 = n[0];
        double n2 = n[2];
        double n3 = n[6];
        double n4 = n[8];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE4(n4,n2,n3,n1,all,r,accuracy);// n1 <-> n4 
        *(result + 1) = VAR4(n4,n2,n3,n1,all,*result);// n1 <-> n4
        *(result + 2) = LRT4(n4,n2,n3,n1,all,*result);// n1 <-> n4
}

//This Headfiel is for group "P1BC2RIL"
double MLE5(double n1,double n2,double n3,double n4,double n,double r0,double accuracy){
        double d1;//[ln(L(r))]'
	double d2;//[ln(L(r))]''
	double r = r0;
	double temp;
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的

	double t1,t2,t3,t4;
	do{
		t1 = 10+2*r;
		t2 = 7+10*r+r*r;
		t3 = 4-2*r;
		t4 = 4*r-r*r;
		d1 = n1*t1/t2 + (n2+n3)*t3/t4 + (-2)*n4/(1-r) - 2*n/(1+2*r);
		d2 = n1*(2*t2-t1*t1)/(t2*t2) + (n2+n3)*((-2)*t4-t3*t3)/(t4*t4) + 2*n4/(pow((1-r),2)) - 4*n/(pow((1+2*r),2));
		temp = d1/d2;
		r -=temp;
                num++;
                if(fabs(temp) < min){min = fabs(temp);minR = r;}
                if(num > 1000){r = minR;break;}

	}while(fabs(temp)>accuracy);
	r = corrcteMLE(r);
	return r;
}

double VAR5(double n1,double n2,double n3,double n4,double n,double r){
	double z;
	double t1,t2,t3,t4;
	t1 = 10+2*r;
	t2 = 7+10*r+r*r;
	t3 = 4-2*r;
	t4 = 4*r-r*r;
	z = n1*(2*t2-t1*t1)/(t2*t2) + (n2+n3)*((-2)*t4-t3*t3)/(t4*t4) + 2*n4/(pow((1-r),2)) - 4*n/(pow((1+2*r),2));
	return (-1)/z;
}

double LRT5(double n1,double n2,double n3,double n4,double n,double R){
        double r,a,b,z;
	r = R;
	a = n1*log(7+10*r+r*r) + (n2+n3)*log(4*r-r*r) + 2*n4*log(1-r) - n*log(1+2*r);
	r = 0.5;
	b = n1*log(7+10*r+r*r) + (n2+n3)*log(4*r-r*r) + 2*n4*log(1-r) - n*log(1+2*r);
	z = (-2)*(b-a);
	return z;
}

void p1bc2ril(float n[],double *result,double r,double accuracy){
        double n1 = n[0];
        double n2 = n[2];
        double n3 = n[6];
        double n4 = n[8];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE5(n1,n2,n3,n4,all,r,accuracy);
        *(result + 1) = VAR5(n1,n2,n3,n4,all,*result);
        *(result + 2) = LRT5(n1,n2,n3,n4,all,*result);
}

//This is for group "P2BC2RIL"
double MLE6(double n1,double n2,double n3,double n4,double n,double r0,double accuracy){
        double d1;//[ln(L(r))]'
	double d2;//[ln(L(r))]''
	double r = r0;
	double temp;
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的
	double t1,t2,t3,t4;
	do{
		t1 = 10+2*r;
		t2 = 7+10*r+r*r;
		t3 = 4-2*r;
		t4 = 4*r-r*r;
		d1 = n1*t1/t2 + (n2+n3)*t3/t4 + (-2)*n4/(1-r) - 2*n/(1+2*r);
		d2 = n1*(2*t2-t1*t1)/(t2*t2) + (n2+n3)*((-2)*t4-t3*t3)/(t4*t4) + 2*n4/(pow((1-r),2)) - 4*n/(pow((1+2*r),2));
		temp = d1/d2;
		r -=temp;
                num++;
                if(fabs(temp) < min){min = fabs(temp);minR = r;}
                if(num > 1000){r = minR;break;}
	}while(fabs(temp)>accuracy);
	r = corrcteMLE(r);
	return r;
}

double VAR6(double n1,double n2,double n3,double n4,double n,double r){
	double z;
	double t1,t2,t3,t4;
	t1 = 10+2*r;
	t2 = 7+10*r+r*r;
	t3 = 4-2*r;
	t4 = 4*r-r*r;
	z = n1*(2*t2-t1*t1)/(t2*t2) + (n2+n3)*((-2)*t4-t3*t3)/(t4*t4) + 2*n4/(pow((1-r),2)) - 4*n/(pow((1+2*r),2));
	return (-1)/z;
}

double LRT6(double n1,double n2,double n3,double n4,double n,double R){
        double r,a,b,z;
	r = R;
	a = n1*log(7+10*r+r*r) + (n2+n3)*log(4*r-r*r) + 2*n4*log(1-r) - n*log(1+2*r);
	r = 0.5;
	b = n1*log(7+10*r+r*r) + (n2+n3)*log(4*r-r*r) + 2*n4*log(1-r) - n*log(1+2*r);
	z = (-2)*(b-a);
	return z;
}

void p2bc2ril(float n[],double *result,double r,double accuracy){
        double n1 = n[8];//n1 <-> n4
        double n2 = n[2];
        double n3 = n[6];
        double n4 = n[0];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE6(n1,n2,n3,n4,all,r,accuracy);
        *(result + 1) = VAR6(n1,n2,n3,n4,all,*result);
        *(result + 2) = LRT6(n1,n2,n3,n4,all,*result);
}

//This is for group "P1BC1DH"
double MLE7(double n1,double n2,double n3,double n4,double n,double r0,double accuracy){
        double b,t,r;
	b = 2*n - 3*n1 - n4;
	t = (sqrt(b*b+8*n*n4)-b) / (2*n);
	r = 1 - sqrt(t);
	r = corrcteMLE(r);
	return r;
}

double VAR7(double n1,double n2,double n3,double n4,double n,double r){
	//t = (1 - r)^2
	//t1 = 2r -2
	//t2 = 2
	double t,t1,t2,i,v;
	t = (1 - r)*(1 - r);
	t1 = 2*r - 2;
	t2 = 2;
	i = n1*(t2*(2+t)-t1*t1)/(pow((2+t),2)) + (n2+n3)*(t2*(t-1)-t1*t1)/(pow((1-t),2)) + n4*(t2*t-t1*t1)/(pow(t,2));
	v = (-1)/i;
	return v;
}

double LRT7(double n1,double n2,double n3,double n4,double n,double R){
        double a,b,t,r,lrt;
	r = 0.5;
	t = pow((1-r),2);	
	a = n1*log(2+t) + (n2+n3)*log(1-t) + n4*log(t);
	r = R;
	t = pow((1-r),2);
	b = n1*log(2+t) + (n2+n3)*log(1-t) + n4*log(t);
	lrt = (-2) * (a - b);
	return lrt;
}

void p1bc1dh(float n[],double *result,double r,double accuracy){
        double n1 = n[0];
        double n2 = n[2];
        double n3 = n[6];
        double n4 = n[8];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE7(n1,n2,n3,n4,all,r,accuracy);
        *(result + 1) = VAR7(n1,n2,n3,n4,all,*result);
        *(result + 2) = LRT7(n1,n2,n3,n4,all,*result);
}

//This Headfiel is for group "P2BC1DH"
double MLE8(double n1,double n2,double n3,double n4,double n,double r0,double accuracy){
        double b,t,r;
	b = 2*n - 3*n1 - n4;
	t = (sqrt(b*b+8*n*n4)-b) / (2*n);
	r = 1 - sqrt(t);
	r = corrcteMLE(r);
	return r;
}

double VAR8(double n1,double n2,double n3,double n4,double n,double r){
	//t = (1 - r)^2
	//t1 = 2r -2
	//t2 = 2
	double t,t1,t2,i,v;
	t = (1 - r)*(1 - r);
	t1 = 2*r - 2;
	t2 = 2;
	i = n1*(t2*(2+t)-t1*t1)/(pow((2+t),2)) + (n2+n3)*(t2*(t-1)-t1*t1)/(pow((1-t),2)) + n4*(t2*t-t1*t1)/(pow(t,2));
	v = (-1)/i;
	return v;
}

double LRT8(double n1,double n2,double n3,double n4,double n,double R){
        double a,b,t,r,lrt;
	r = 0.5;
	t = pow((1-r),2);	
	a = n1*log(2+t) + (n2+n3)*log(1-t) + n4*log(t);
	r = R;
	t = pow((1-r),2);
	b = n1*log(2+t) + (n2+n3)*log(1-t) + n4*log(t);
	lrt = (-2) * (a - b);
	return lrt;
}

void p2bc1dh(float n[],double *result,double r,double accuracy){
        double n4 = n[0];
        double n2 = n[2];
        double n3 = n[6];
        double n1 = n[8];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE8(n1,n2,n3,n4,all,r,accuracy);
        *(result + 1) = VAR8(n1,n2,n3,n4,all,*result);
        *(result + 2) = LRT8(n1,n2,n3,n4,all,*result);
}

//This is for group "P1BC2DH"
double MLE9(double n1,double n2,double n3,double n4,double n,double r0,double accuracy){
        double b,t,r;
	b = 6*n - 7*n1 - n4;
	t = (sqrt(b*b + 24*n*n4)-b) / (2*n);
	r = 1 - pow(t,1.0/3);
	r = corrcteMLE(r);
	return r;
}

double VAR9(double n1,double n2,double n3,double n4,double n,double r){
	//t = (1 - r)^3
	//t1 = -3(1-r)^2
	//t2 = 6-6r
	double t,t1,t2,i,v;
	t = (1 - r)*(1 - r)*(1 - r);
	t1 = (-3) *(1 - r)*(1 - r);
	t2 = 6 - 6*r;
	i = n1*(t2*(6+t)-t1*t1)/(pow((6+t),2)) + (n2+n3)*(t2*(t-1)-t1*t1)/(pow((1-t),2)) + n4*(t2*t-t1*t1)/(pow(t,2));
	v = (-1)/i;
	return v;
}

double LRT9(double n1,double n2,double n3,double n4,double n,double R){
        double a,b,t,r,lrt;
	r = 0.5;
	t = pow((1-r),3);	
	a = n1*log(6+t) + (n2+n3)*log(1-t) + n4*log(t);
	r = R;
	t = pow((1-r),3);
	b = n1*log(6+t) + (n2+n3)*log(1-t) + n4*log(t);
	lrt = (-2) * (a - b);
	return lrt;
}

void p1bc2dh(float n[],double *result,double r,double accuracy){
        double n1 = n[0];
        double n2 = n[2];
        double n3 = n[6];
        double n4 = n[8];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE9(n1,n2,n3,n4,all,r,accuracy);
        *(result + 1) = VAR9(n1,n2,n3,n4,all,*result);
        *(result + 2) = LRT9(n1,n2,n3,n4,all,*result);
}

//This is for group "P2BC2DH"
double MLE10(double n1,double n2,double n3,double n4,double n,double r0,double accuracy){
        double b,t,r;
	b = 6*n - 7*n1 - n4;
	t = (sqrt(b*b + 24*n*n4)-b) / (2*n);
	r = 1 - pow(t,1.0/3);
	r = corrcteMLE(r);
	return r;
}

double VAR10(double n1,double n2,double n3,double n4,double n,double r){
	//t = (1 - r)^3
	//t1 = -3(1-r)^2
	//t2 = 6-6r
	double t,t1,t2,i,v;
	t = (1 - r)*(1 - r)*(1 - r);
	t1 = (-3) *(1 - r)*(1 - r);
	t2 = 6 - 6*r;
	i = n1*(t2*(6+t)-t1*t1)/(pow((6+t),2)) + (n2+n3)*(t2*(t-1)-t1*t1)/(pow((1-t),2)) + n4*(t2*t-t1*t1)/(pow(t,2));
	v = (-1)/i;
	return v;
}

double LRT10(double n1,double n2,double n3,double n4,double n,double R){
        double a,b,t,r,lrt;
	r = 0.5;
	t = pow((1-r),3);	
	a = n1*log(6+t) + (n2+n3)*log(1-t) + n4*log(t);
	r = R;
	t = pow((1-r),3);
	b = n1*log(6+t) + (n2+n3)*log(1-t) + n4*log(t);
	lrt = (-2) * (a - b);
	return lrt;
}

void p2bc2dh(float n[],double *result,double r,double accuracy){
        double n1 = n[8];
        double n2 = n[2];
        double n3 = n[6];
        double n4 = n[0];
        double all  = n1 + n2 + n3 + n4;

        *result = MLE10(n1,n2,n3,n4,all,r,accuracy);
        *(result + 1) = VAR10(n1,n2,n3,n4,all,*result);
        *(result + 2) = LRT10(n1,n2,n3,n4,all,*result);
}

//This is for group "P1BC1F1.h"
double MLE11(double n2,double n4,double n){
        double r;
	r = (n2 + n4) / n;
	r = corrcteMLE(r);
	return r;
}

double VAR11(double n1,double n2,double n4,double n5,double r){
	double t1,t2,i,v;
	t1 = (1 - r)*(1 -r);
	t2 = r * r;
	i = (n1 + n5) / t1 + (n2 + n4) / t2;
	v = (-1)/i;
	return v;
}

double LRT11(double n1,double n2,double n4,double n5,double n,double r){
        double lrt;
	lrt = (-2) * (n * log(0.5) - (n1+n5) * log(1 - r) - (n2+n4)*log(r));
	return lrt;
}

void p1bc1f1(float n[],double *result,double r,double accuracy){
        double n1 = n[0];
        double n2 = n[1];
        double n4 = n[3];
        double n5 = n[4];
        double all  = n1 + n2 + n4 + n5;

        *result = MLE11(n2,n4,all);
        *(result + 1) = VAR11(n1,n2,n4,n5,*result);
        *(result + 2) = LRT11(n1,n2,n4,n5,all,*result);
}

//This is for group "P2BC1F1.h"
double MLE12(double n6,double n8,double n){
        double r;
	r = (n6 + n8) / n;
	r = corrcteMLE(r);
	return r;
}

double VAR12(double n5,double n6,double n8,double n9,double n,double r){
	double t1,t2,i,v;
	t1 = (1 - r)*(1 -r);
	t2 = r * r;
	i = (n5 + n9) / t1 + (n6 + n8) / t2;
	v = (-1)/i;
	return v;
}

double LRT12(double n5,double n6,double n8,double n9,double n,double r){
        double lrt;
	lrt = (-2) * (n * log(0.5) - (n5+n9) * log(1 - r) - (n6+n8)*log(r));
	return lrt;
}

void p2bc1f1(float n[],double *result,double r,double accuracy){
        double n5 = n[4];
        double n6 = n[5];
        double n8 = n[7];
        double n9 = n[8];
        double all  = n5 + n6 + n8 + n9;

        *result = MLE12(n6,n8,all);
        *(result + 1) = VAR12(n5,n6,n8,n9,all,*result);
        *(result + 2) = LRT12(n5,n6,n8,n9,all,*result);
}

//This is for group "F2.h"
double MLE13(double n5,double N1,double N2,double r,double accuracy){
	double R = 0.01;
        double d1;//[ln(L(r))]'
        double d2;//[ln(L(r))]''
        double temp;
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的

	do{	
		double t1 = 2*r*r - 2*r + 1;
		double t2 = 4*r - 2;
		d1 = (-N1)/(1-r) + N2/r + n5*t2/t1;
		d2 = (-N1)/(pow((1- r),2)) - N2/pow(r,2) + n5*(4*t1 - t2*t2)/pow(t1,2);
		temp = d1/d2;
		if(isnan(temp) || (num > 20 && R > accuracy)){//temp不能为nan值
			temp = 1;
			R -= (accuracy * 0.1);
			r = R ;
			num = 0;
			min = 1;
		}else{
			r -=temp;
			if(r < 0){temp = 1;}//不能马上输出
			num++;
			if(fabs(temp) < min && r > 0){min = fabs(temp);minR = r;}
		}
                if(num > 21){r = minR;break;}//如果备选的初始值都试过了，则不再尝试
	}while(fabs(temp)>accuracy);
	//r = corrcteMLE(r);
	return r;
}

double VAR13(double n5,double N1,double N2,double r){
	double t1,i,v;
	t1 = 2*r*r + 2*r + 1;
	i = N1/(pow((1- r),2)) + N2/pow(r,2) + 8*n5*r*(r+1)/(t1*t1);
	v = (-1)/i;
	return v;
}

double LRT13(double n5,double N1,double N2,double n,double r){
        double lrt;
	lrt = (-2) * (2*n-n5)*log(0.5) - N1*log(1-r) - N2*log(r) - n5*log(2*r*r + 2*r + 1);
	return lrt;
}


//this function return "MLE" "variance" "LRT" "LOD"  of f2
void f2(float n[],double *result,double r,double accuracy){
        double n1 = n[0];
        double n2 = n[1];
        double n3 = n[2];
        double n4 = n[3];
        double n5 = n[4];
        double n6 = n[5];
        double n7 = n[6];
        double n8 = n[7];
        double n9 = n[8];
        double all  = n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8 +n9;
	double N1 = 2 * n1 + n2 + n4 + n6 + n8 + 2 * n9;
	double N2 = n2 + 2 * n3 + n4 + n6 + 2 * n7 + n8;
	

        *result = MLE13(n5,N1,N2,r,accuracy);
        *(result + 1) = VAR13(n5,N1,N2,*result);
        *(result + 2) = LRT13(n5,N1,N2,all,*result);
}

//This is for group "F3.h"
void f3g(float n[],double *result,double R,double accuracy){
        double n1 = n[0];
        double n2 = n[1];
        double n3 = n[2];
        double n4 = n[3];
        double n5 = n[4];
        double n6 = n[5];
        double n7 = n[6];
        double n8 = n[7];
	double n9 = n[8];

	double f1,f1d1,f1d2;//f1, f1d1=(f1)' f1d2=(f1)''
	double f2,f2d1,f2d2;
	double f3,f3d1,f3d2;
	double f5,f5d1,f5d2;

	double r = R;//initialize r
	double r4,r3,r2;//r^4,r^3,r^2
	double d1,d2,temp;//d1 = lnL(r)';d2 = lnL(r)'' ;temp = d1/d2;
	float flag = 0;//temp less than accuracy or not
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的
	do{
		r2 = r * r;
		r3 = r2 * r;
		r4 = r3 * r;
	
		f1 = (2*r4 - 4*r3 + 6*r2 - 6*r + 3)/8;
		f1d1 = (8*r3 - 12*r2 + 12*r - 6 )/8;
		f1d2 = (24*r2 - 24*r + 12)/8;

		f2 = ((-1) * r4 + 2*r3 - 2*r2 + r)/2;
		f2d1 = ((-4)*r3 + 6*r2 -4*r + 1)/2;
		f2d2 = ((-12)*r2 + 12*r - 4)/2;

		f3 = (r4 - 2*r3 + r2 + r)/4;
		f3d1 = (4*r3 - 6*r2 +2*r + 1)/4;
		f3d2 = (12*r2 - 12*r + 2)/4;

		f5 = (4*r4 - 8*r3 + 8*r2 - 4*r + 1)/4;
		f5d1 = (16*r3 - 24*r2 + 16*r -4)/4;
		f5d2 = (48*r2 + 48*r + 16)/4;
	
		d1 = (n1 + n9) * f1d1 / f1 + (n2 + n4 + n6 + n8) * f2d1 / f2 + (n3 + n7) * f3d1 / f3 + n5 * f5d1 / f5;
		d2 = (n1 + n9) * (f1d2 * f1 - f1d1 * f1d1) / (f1 * f1) + (n2 + n4 + n6 + n8) * (f2d2 * f2 - f2d1 * f2d1) / (f2 * f2) + (n3 + n7) * (f3d2 * f3 - f3d1 * f3d1) / (f3 * f3) + n5 * (f5d2 * f5 - f5d1 * f5d1) / (f5 * f5);
		if(flag == 1){break;}
		
		temp = d1 / d2;
		r = r - temp;
                num++;
                if(fabs(temp) < min){min = fabs(temp);minR = r;}
                if(num > 1000){r = minR;flag = 1;}
		if( fabs(temp) <= accuracy){flag = 1;}//if "temp" less than "accuracy", "r" is finish, and calculated "lnL(r)''" last times; 

	}while(1);
	r = corrcteMLE(r);
        *result = r;//MLE
        *(result + 1) = -1.0 / d2;//VAR

	double x,y;
	y = (n1 + n9) * log(f1) + (n2 + n4 + n6 + n8) * log(f2) + (n3 + n7) * log(f3) + n5 * log(f5);
	r = 0.5;
	r2 = r * r;
	r3 = r2 * r;
	r4 = r3 * r;
	f1 = (2*r4 - 4*r3 + 6*r2 - 6*r + 3)/8;
	f2 = ((-1) * r4 + 2*r3 - 2*r2 + r)/2;
	f3 = (r4 - 2*r3 + r2 + r)/4;
	f5 = (4*r4 - 8*r3 + 8*r2 - 4*r + 1)/4;
	x = (n1 + n9) * log(f1) + (n2 + n4 + n6 + n8) * log(f2) + (n3 + n7) * log(f3) + n5 * log(f5);

	
        *(result + 2) = -2 * (x - y);//LRT
}

//This is for group "P1BC2F1"
double MLE15(double n1,double n5,double n){
        double b,t,r;
	b = 2*n - 3*n1 - n5;
	t = (sqrt(b*b + 8*n*n5)-b) / (2*n);
	r = 1 - sqrt(t);
	r = corrcteMLE(r);
	return r;
}

double VAR15(double n1,double n2,double n4,double n5,double n,double r){
	//t = (1 - r)^2
	//t1 = -2(1-r)
	//t2 = 2
	double t,t1,t2,i,v;
	t = (1 - r)*(1 - r);
	t1 = (-2) *(1 - r);
	t2 = 2;
	i = n1*(t2*(2+t)-t1*t1)/(pow((2+t),2)) - (n2+n4)*(t2*(1-t)-t1*t1)/(pow((1-t),2)) + n5*(t2*t-t1*t1)/(pow(t,2));
	v = (-1)/i;
	return v;
}

double LRT15(double n1,double n2,double n4,double n5,double R){
        double a,b,t,r,lrt;
	r = 0.5;
	t = pow((1-r),2);	
	a = n1*log(2+t) + (n2+n4)*log(1-t) + n5*log(t);
	r = R;
	t = pow((1-r),2);
	b = n1*log(2+t) + (n2+n4)*log(1-t) + n5*log(t);
	lrt = (-2) * (a - b);
	return lrt;
}

void p1bc2f1(float n[],double *result,double r,double accuracy){
        double n1 = n[0];
        double n2 = n[1];
        double n4 = n[3];
	double n5 = n[4];
        double all  = n1 + n2 + n4 + n5;

        *result = MLE15(n1,n5,all);
        *(result + 1) = VAR15(n1,n2,n4,n5,all,*result);
        *(result + 2) = LRT15(n1,n2,n4,n5,*result);
}

//This is for group "P2BC2F1.h"
double MLE16(double n5,double n9,double n){
        double b,t,r;
        b = 2*n - 3*n9 - n5;
        t = (sqrt(b*b + 8*n*n5)-b) / (2*n);
        r = 1 - sqrt(t);
	r = corrcteMLE(r);
        return r;
}

double VAR16(double n5,double n6,double n8,double n9,double n,double r){
        //t = (1 - r)^2
        //t1 = -2(1-r)
        //t2 = 2
        double t,t1,t2,i,v;
        t = (1 - r)*(1 - r);
        t1 = (-2) *(1 - r);
        t2 = 2;
        i = n9*(t2*(2+t)-t1*t1)/(pow((2+t),2)) - (n6+n8)*(t2*(1-t)-t1*t1)/(pow((1-t),2)) + n5*(t2*t-t1*t1)/(pow(t,2));
        v = (-1)/i;
        return v;
}

double LRT16(double n5,double n6,double n8,double n9,double R){
        double a,b,t,r,lrt;
        r = 0.5;
        t = pow((1-r),2);
        a = n9*log(2+t) + (n6+n8)*log(1-t) + n5*log(t);
        r = R;
        t = pow((1-r),2);
        b = n9*log(2+t) + (n6+n8)*log(1-t) + n5*log(t);
        lrt = (-2) * (a - b);
        return lrt;
}

void p2bc2f1(float n[],double *result,double r,double accuracy){
        double n5 = n[4];
        double n6 = n[5];
        double n8 = n[7];
        double n9 = n[8];
        double all  = n5 + n6 + n8 + n9;

        *result = MLE16(n5,n9,all);
        *(result + 1) = VAR16(n5,n6,n8,n9,all,*result);
        *(result + 2) = LRT16(n5,n6,n8,n9,*result);
}


//This Headfiel is for group "P1BC1F2"
void p1bc1f2(float n[],double *result,double R,double accuracy){
        double n1 = n[0];
        double n2 = n[1];
        double n3 = n[2];
        double n4 = n[3];
        double n5 = n[4];
        double n6 = n[5];
        double n7 = n[6];
        double n8 = n[7];
	double n9 = n[8];

	double f1,f1d1,f1d2;//f1, f1d1=(f1)' f1d2=(f1)''
	double f2,f2d1,f2d2;
	double f3,f3d1,f3d2;
	double f4,f4d1,f4d2;
	double f5,f5d1,f5d2;
	double f6,f6d1,f6d2;
	double f9,f9d1,f9d2;

	double r = R;//initialize r
	double r3,r2;//r^3,r^2
	double d1,d2,temp;//d1 = lnL(r)';d2 = lnL(r)'' ;temp = d1/d2;
	float flag = 0;//temp less than accuracy or not
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的
	do{
		r2 = r * r;
		r3 = r2 * r;
	
		f1 = ( (-1)*r3 + 3*r2 - 5*r + 3)/8;
		f1d1 = ( (-3)*r2 + 6*r - 5 )/8;
		f1d2 = ((-6)*r + 6)/8;

		f2 = (r2 - r + 1)/4;
		f2d1 = (2*r -1)/4;
		f2d2 = 0.5;
		f3 = ( (-1)*r3 + r2 + r)/8;
		f3d1 = ((-3)*r2 +2*r + 1)/8;
		f3d2 = ((-6)*r + 2)/8;

		f4 = (r3 - 2*r2 + 2*r)/4;
		f4d1 = (3*r2 - 4*r +2)/4;
		f4d2 = (6*r - 4)/4;

		f5 = ((-2)*r3 + 4*r2 - 3*r + 1)/4;
		f5d1 = ((-6)*r2 + 8*r -3)/4;
		f5d2 = ((-12)*r + 8)/4;

		f6 = (r3 - 2*r2 + r)/4;
		f6d1 = (3*r2 - 4*r +1)/4;
		f6d2 = (6*r -4)/4;

		f9 = ((-1)*r3 + 3*r2 - 3*r + 1)/8;
		f9d1 = ((-3)*r2 + 6*r -3)/8;
		f9d2 = ((-6)*r + 6)/8;
		
		d1 = n1 * f1d1 / f1 + n2 * f2d1 / f2 + (n3+n7) * f3d1 / f3 + n4 * f4d1 / f4 + n5* f5d1 / f5 + (n6+n8) * f6d1 / f6 + n9 * f9d1 / f9;
		d2 = n1 * (f1d2 * f1 - f1d1 * f1d1) / (f1 * f1) + n2 * (f2d2 * f2 - f2d1 * f2d1) / (f2 * f2) + (n3 + n7) * (f3d2 * f3 - f3d1 * f3d1) / (f3 * f3) + n4 * (f4d2 * f4 - f4d1 * f4d1) / (f4 * f4) + n5 * (f5d2 * f5 - f5d1 * f5d1) / (f5 * f5) + (n6+n8) * (f6d2 * f6 - f6d1 * f6d1) / (f6 * f6) + n9 * (f9d2 * f9 - f9d1 * f9d1) / (f9 * f9);
		if(flag == 1){break;}
		
		temp = d1 / d2;
		r = r - temp;
                num++;
                if(fabs(temp) < min){min = fabs(temp);minR = r;}
                if(num > 1000){r = minR;flag = 1;}

		if( fabs(temp) <= accuracy){flag = 1;}//if "temp" less than "accuracy", "r" is finish, and calculated "lnL(r)'',f1,f2,f3..." last times; 
 
	}while(1);
	//prfloatf("finish\n");
	r = corrcteMLE(r);
        *result = r;//MLE
        *(result + 1) = -1.0 / d2;//VAR

	double x,y;
	y = n1 * log(f1) + n2 * log(f2) + (n3 + n7) * log(f3) + n4 * log(f4) + n5 * log(f5) + (n6+n8) * log(f6) + n9 * log(f9);
	r = 0.5;
	r2 = r * r;
	r3 = r2 * r;
	f1 = ( (-1)*r3 + 3*r2 - 5*r + 3)/8;
	f2 = (r2 - r + 1)/4;
	f3 = ( (-1)*r3 + r2 + r)/8;
	f4 = (r3 - 2*r2 + 2*r)/4;
	f5 = ((-2)*r3 + 4*r2 - 3*r + 1)/4;
	f6 = (r3 - 2*r2 + r)/4;
	f9 = ((-1)*r3 + 3*r2 - 3*r + 1)/8;
	x = n1 * log(f1) + n2 * log(f2) + (n3 + n7) * log(f3) + n4 * log(f4) + n5 * log(f5) + (n6+n8) * log(f6) + n9 * log(f9);

        *(result + 2) = -2 * (x - y);//LRT
}

//This Headfiel is for group "P2BC1F2"
void p2bc1f2(float n[],double *result,double R,double accuracy){
        double n1 = n[0];
        double n2 = n[1];
        double n3 = n[2];
        double n4 = n[3];
        double n5 = n[4];
        double n6 = n[5];
        double n7 = n[6];
        double n8 = n[7];
	double n9 = n[8];

	double f1,f1d1,f1d2;//f1, f1d1=(f1)' f1d2=(f1)''
	double f2,f2d1,f2d2;
	double f3,f3d1,f3d2;
	double f8,f8d1,f8d2;
	double f5,f5d1,f5d2;
	double f6,f6d1,f6d2;
	double f9,f9d1,f9d2;

	double r = R;//initialize r
	double r3,r2;//r^3,r^2
	double d1,d2,temp;//d1 = lnL(r)';d2 = lnL(r)'' ;temp = d1/d2;
	float flag = 0;//temp less than accuracy or not
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的
	do{
		r2 = r * r;
		r3 = r2 * r;
	
		f9 = ((-1)*r3 + 3*r2 - 5*r + 3)/8;
		f9d1 = ( (-3)*r2 + 6*r - 5 )/8;
		f9d2 = ((-6)*r + 6)/8;

		f8 = (r2 - r + 1)/4;
		f8d1 = (2*r -1)/4;
		f8d2 = 0.5;

		f3 = ( (-1)*r3 + r2 + r)/8;
		f3d1 = ((-3)*r2 +2*r + 1)/8;
		f3d2 = ((-6)*r + 2)/8;

		f6 = (r3 - 2*r2 + 2*r)/4;
		f6d1 = (3*r2 - 4*r +2)/4;
		f6d2 = (6*r - 4)/4;

		f5 = ((-2)*r3 + 4*r2 - 3*r + 1)/4;
		f5d1 = ((-6)*r2 + 8*r -3)/4;
		f5d2 = ((-12)*r + 8)/4;

		f2 = (r3 - 2*r2 + r)/4;
		f2d1 = (3*r2 - 4*r +1)/4;
		f2d2 = (6*r -4)/4;

		f1 = ((-1)*r3 + 3*r2 - 3*r + 1)/8;
		f1d1 = ((-3)*r2 + 6*r -3)/8;
		f1d2 = ((-6)*r + 6)/8;
		
		d1 = n1 * f1d1 / f1 + (n2+n4) * f2d1 / f2 + (n3+n7) * f3d1 / f3 + n5* f5d1 / f5 + n6 * f6d1 / f6 + n8 * f8d1 / f8 + n9 * f9d1 / f9;
		d2 = n1 * (f1d2 * f1 - f1d1 * f1d1) / (f1 * f1) + (n2+n4) * (f2d2 * f2 - f2d1 * f2d1) / (f2 * f2) + (n3 + n7) * (f3d2 * f3 - f3d1 * f3d1) / (f3 * f3) + n5 * (f5d2 * f5 - f5d1 * f5d1) / (f5 * f5) + n6 * (f6d2 * f6 - f6d1 * f6d1) / (f6 * f6) +  n8 * (f8d2 * f8 - f8d1 * f8d1) / (f8 * f8) + n9 * (f9d2 * f9 - f9d1 * f9d1) / (f9 * f9);
		if(flag == 1){break;}
		
		temp = d1 / d2;
		r = r - temp;
                num++;
                if(fabs(temp) < min){min = fabs(temp);minR = r;}
                if(num > 1000){r = minR;flag = 1;}

		if( fabs(temp) <= accuracy){flag = 1;}//if "temp" less than "accuracy", "r" is finish, and calculated "lnL(r)'',f1,f2,f3..." last times; 
 
	}while(1);
	//prfloatf("finish\n");
	r = corrcteMLE(r);
        *result = r;//MLE
        *(result + 1) = -1.0 / d2;//VAR

	double x,y;
	y = n1 * log(f1) + (n2+n4) * log(f2) + (n3 + n7) * log(f3) + n5 * log(f5) + n6 * log(f6) + n8 * log(f8) + n9 * log(f9);
	r = 0.5;
	r2 = r * r;
	r3 = r2 * r;
	f1 = ((-1)*r3 + 3*r2 - 3*r + 1)/8;
	f2 = (r3 - 2*r2 + r)/4;
	f3 = ( (-1)*r3 + r2 + r)/8;
	f5 = ((-2)*r3 + 4*r2 - 3*r + 1)/4;
	f6 = (r3 - 2*r2 + 2*r)/4;
	f8 = (r2 - r + 1)/4;
	f9 = ((-1)*r3 + 3*r2 - 5*r + 3)/8;
	x = n1 * log(f1) + (n2+n4) * log(f2) + (n3 + n7) * log(f3) + n5 * log(f5) + n6 * log(f6) + n8 * log(f8) + n9 * log(f9);
        *(result + 2) = -2 * (x - y);//LRT
}
//This is for group "P1BC2F2"
void p1bc2f2(float n[],double *result,double R,double accuracy){
        double n1 = n[0];
        double n2 = n[1];
        double n3 = n[2];
        double n4 = n[3];
        double n5 = n[4];
        double n6 = n[5];
        double n7 = n[6];
        double n8 = n[7];
	double n9 = n[8];

	double f1,f1d1,f1d2;//f1, f1d1=(f1)' f1d2=(f1)''
	double f2,f2d1,f2d2;
	double f3,f3d1,f3d2;
	double f5,f5d1,f5d2;
	double f6,f6d1,f6d2;
	double f9,f9d1,f9d2;

	double r = R;//initialize r
	double r4,r3,r2;//r^3,r^2
	double d1,d2,temp;//d1 = lnL(r)';d2 = lnL(r)'' ;temp = d1/d2;
	float flag = 0;//temp less than accuracy or not
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的
	do{
		r2 = r * r;
		r3 = r2 * r;
		r4 = r3 * r;
		f1 = (r4 - 4*r3 + 8*r2 - 8*r + 13)/16;
		f1d1 = ( 4*r3 - 12*r2 + 16*r - 8 )/16;
		f1d2 = (12*r2 - 24*r + 16)/16;

		f2 = ((-1)*r4 + 3*r3 - 4*r2 + 3*r)/8;
		f2d1 = ((-4)*r3 + 9*r2 - 8*r + 3)/8;
		f2d2 = ((-12)*r2 + 18*r - 8)/8;

		f3 = (r4 - 2*r3 + 2*r)/16;
		f3d1 = (4*r3 - 6*r2 + 2)/16;
		f3d2 = (12*r2 - 12*r)/16;

		f5 = (2*r4 - 6*r3 + 7*r2 - 4*r + 1)/8;
		f5d1 = (8*r3 - 18*r2 + 14*r - 4)/8;
		f5d2 = (24*r2 - 36*r + 14)/8;

		f6 = ((-1)*r4 + 3*r3 - 3*r2 + r)/8;
		f6d1 = ((-4)*r3 + 9*r2 - 6*r +1)/8;
		f6d2 = ((-12)*r2 + 18*r - 6)/8;

		f9 = (r4 - 4*r3 + 6*r2 - 4*r + 1)/16;
		f9d1 = (4*r3 - 12*r2 + 12*r - 4)/8;
		f9d2 = (12*r2 - 24*r + 12)/8;
		
		d1 = n1 * f1d1 / f1 + (n2+n4) * f2d1 / f2 + (n3+n7) * f3d1 / f3 + n5* f5d1 / f5 + (n6+n8) * f6d1 / f6 + n9 * f9d1 / f9;
		d2 = n1 * (f1d2 * f1 - f1d1 * f1d1) / (f1 * f1) + (n2+n4) * (f2d2 * f2 - f2d1 * f2d1) / (f2 * f2) + (n3 + n7) * (f3d2 * f3 - f3d1 * f3d1) / (f3 * f3) + n5 * (f5d2 * f5 - f5d1 * f5d1) / (f5 * f5) + (n6+n8) * (f6d2 * f6 - f6d1 * f6d1) / (f6 * f6) + n9 * (f9d2 * f9 - f9d1 * f9d1) / (f9 * f9);
		if(flag == 1){break;}
		
		temp = d1 / d2;
		r = r - temp;
                num++;
                if(fabs(temp) < min){min = fabs(temp);minR = r;}
                if(num > 1000){r = minR;flag = 1;}

		if( fabs(temp) <= accuracy){flag = 1;}//if "temp" less than "accuracy", "r" is finish, and calculated "lnL(r)'',f1,f2,f3..." last times; 
 
	}while(1);
	//prfloatf("finish\n");
	r = corrcteMLE(r);
        *result = r;//MLE
        *(result + 1) = -1.0 / d2;//VAR

	double x,y;
	y = n1 * log(f1) + (n2+n4) * log(f2) + (n3 + n7) * log(f3) + n5 * log(f5) + (n6+n8) * log(f6) + n9 * log(f9);
	r = 0.5;
	r2 = r * r;
	r3 = r2 * r;
	r4 = r3 * r;
	f1 = (r4 - 4*r3 + 8*r2 - 8*r + 13)/16;
	f2 = ((-1)*r4 + 3*r3 - 4*r2 + 3*r)/8;
	f3 = (r4 - 2*r3 + 2*r)/16;
	f5 = (2*r4 - 6*r3 + 7*r2 - 4*r + 1)/8;
	f6 = ((-1)*r4 + 3*r3 - 3*r2 + r)/8;
	f9 = (r4 - 4*r3 + 6*r2 - 4*r + 1)/16;
	x = n1 * log(f1) + (n2+n4) * log(f2) + (n3 + n7) * log(f3) + n5 * log(f5) + (n6+n8) * log(f6) + n9 * log(f9);

        *(result + 2) = -2 * (x - y);//LRT
}

//This is for group "P2BC2F2"
void p2bc2f2(float n[],double *result,double R,double accuracy){
        double n1 = n[0];
        double n2 = n[1];
        double n3 = n[2];
        double n4 = n[3];
        double n5 = n[4];
        double n6 = n[5];
        double n7 = n[6];
        double n8 = n[7];
	double n9 = n[8];

	double f1,f1d1,f1d2;//f1, f1d1=(f1)' f1d2=(f1)''
	double f2,f2d1,f2d2;
	double f3,f3d1,f3d2;
	double f5,f5d1,f5d2;
	double f6,f6d1,f6d2;
	double f9,f9d1,f9d2;

	double r = R;//initialize r
	double r4,r3,r2;//r^3,r^2
	double d1,d2,temp;//d1 = lnL(r)';d2 = lnL(r)'' ;temp = d1/d2;
	float flag = 0;//temp less than accuracy or not
        double min = 1;
        double minR;
        int num = 0;//循环次数 如果循环次数大于1000，则不再循环，选取1000次中最小的
	do{
		r2 = r * r;
		r3 = r2 * r;
		r4 = r3 * r;
		f9 = (r4 - 4*r3 + 8*r2 - 8*r + 13)/16;
		f9d1 = ( 4*r3 - 12*r2 + 16*r - 8 )/16;
		f9d2 = (12*r2 - 24*r + 16)/16;

		f6 = ((-1)*r4 + 3*r3 - 4*r2 + 3*r)/8;
		f6d1 = ((-4)*r3 + 9*r2 - 8*r + 3)/8;
		f6d2 = ((-12)*r2 + 18*r - 8)/8;

		f3 = (r4 - 2*r3 + 2*r)/16;
		f3d1 = (4*r3 - 6*r2 + 2)/16;
		f3d2 = (12*r2 - 12*r)/16;

		f5 = (2*r4 - 6*r3 + 7*r2 - 4*r + 1)/8;
		f5d1 = (8*r3 - 18*r2 + 14*r - 4)/8;
		f5d2 = (24*r2 - 36*r + 14)/8;

		f2 = ((-1)*r4 + 3*r3 - 3*r2 + r)/8;
		f2d1 = ((-4)*r3 + 9*r2 - 6*r +1)/8;
		f2d2 = ((-12)*r2 + 18*r - 6)/8;

		f1 = (r4 - 4*r3 + 6*r2 - 4*r + 1)/16;
		f1d1 = (4*r3 - 12*r2 + 12*r - 4)/8;
		f1d2 = (12*r2 - 24*r + 12)/8;
		
		d1 = n1 * f1d1 / f1 + (n2+n4) * f2d1 / f2 + (n3+n7) * f3d1 / f3 + n5* f5d1 / f5 + (n6+n8) * f6d1 / f6 + n9 * f9d1 / f9;
		d2 = n1 * (f1d2 * f1 - f1d1 * f1d1) / (f1 * f1) + (n2+n4) * (f2d2 * f2 - f2d1 * f2d1) / (f2 * f2) + (n3 + n7) * (f3d2 * f3 - f3d1 * f3d1) / (f3 * f3) + n5 * (f5d2 * f5 - f5d1 * f5d1) / (f5 * f5) + (n6+n8) * (f6d2 * f6 - f6d1 * f6d1) / (f6 * f6) + n9 * (f9d2 * f9 - f9d1 * f9d1) / (f9 * f9);
		if(flag == 1){break;}
		
		temp = d1 / d2;
		r = r - temp;
                num++;
                if(fabs(temp) < min){min = fabs(temp);minR = r;}
                if(num > 1000){r = minR;flag = 1;}

		if( fabs(temp) <= accuracy){flag = 1;}//if "temp" less than "accuracy", "r" is finish, and calculated "lnL(r)'',f1,f2,f3..." last times; 
 
	}while(1);
	//prfloatf("finish\n");
	r = corrcteMLE(r);
        *result = r;//MLE
        *(result + 1) = -1.0 / d2;//VAR

	double x,y;
	y = n1 * log(f1) + (n2+n4) * log(f2) + (n3 + n7) * log(f3) + n5 * log(f5) + (n6+n8) * log(f6) + n9 * log(f9);
	r = 0.5;
	r2 = r * r;
	r3 = r2 * r;
	r4 = r3 * r;
	f9 = (r4 - 4*r3 + 8*r2 - 8*r + 13)/16;
	f6 = ((-1)*r4 + 3*r3 - 4*r2 + 3*r)/8;
	f3 = (r4 - 2*r3 + 2*r)/16;
	f5 = (2*r4 - 6*r3 + 7*r2 - 4*r + 1)/8;
	f2 = ((-1)*r4 + 3*r3 - 3*r2 + r)/8;
	f1 = (r4 - 4*r3 + 6*r2 - 4*r + 1)/16;
	x = n1 * log(f1) + (n2+n4) * log(f2) + (n3 + n7) * log(f3) + n5 * log(f5) + (n6+n8) * log(f6) + n9 * log(f9);

        *(result + 2) = (-2) * (x - y);//LRT
}


//计算LOD
double callLOD(double lrt,double mle){
	if(mle <= 0){
		return 30;
	}else{
		return lrt/(2*log(10));
	}
}


//以下是用于CP群体的

void get_mi_matrix(int phenotype,char mark,double (*mi)[1]){

	//double (*mi)[1] = (double (*)[1])mi_point;
	if(phenotype == 1){//abxcd
		if(mark == '1'){ *(*(mi+0)+0) = 1;*(*(mi+1)+0) = 0;*(*(mi+2)+0) = 0;*(*(mi+3)+0) = 0;}//ac
		else if(mark == '2'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 1;*(*(mi+2)+0) = 0;*(*(mi+3)+0) = 0;}//ad
		else if(mark == '3'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 0;*(*(mi+2)+0) = 1;*(*(mi+3)+0) = 0;}//bc
		else if(mark == '4'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 0;*(*(mi+2)+0) = 0;*(*(mi+3)+0) = 1;}//bd
		else{
			//fprintf(stderr,"Error: abxcd phenotype mark's type must be(ac,ad,bc,bd).\n");
		}
	}else if(phenotype == 2){//abxac
		if(mark == '5'){ *(*(mi+0)+0) = 1;*(*(mi+1)+0) = 0;*(*(mi+2)+0) = 0;*(*(mi+3)+0) = 0;}//aa
                else if(mark == '1'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 1;*(*(mi+2)+0) = 0;*(*(mi+3)+0) = 0;}//ac
                else if(mark == '6'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 0;*(*(mi+2)+0) = 1;*(*(mi+3)+0) = 0;}//ab
                else if(mark == '3'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 0;*(*(mi+2)+0) = 0;*(*(mi+3)+0) = 1;}//bc
                else{
			//fprintf(stderr,"Error: abxac phenotype mark's type must be(aa,ac,ab,bc).\n ");
		}
	}else if(phenotype == 3){//abxab
		if(mark == '5'){ *(*(mi+0)+0) = 1;*(*(mi+1)+0) = 0;*(*(mi+2)+0) = 0;}//aa
                else if(mark == '6'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 1;*(*(mi+2)+0) = 0;}//ab
                else if(mark == '7'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 0;*(*(mi+2)+0) = 1;}//bb
                else{
			//fprintf(stderr,"Error: abxab phenotype mark's type must be(aa,ab,bb).\n ");
		}
	}else if(phenotype == 4){//abxaa
		if(mark == '5'){ *(*(mi+0)+0) = 1;*(*(mi+1)+0) = 0;}//aa
                else if(mark == '6'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 1;}//ab
                else{
			//fprintf(stderr,"Error: abxaa phenotype mark's type must be(aa,ab).\n ");
		}
	}else if(phenotype == 5){//aaxab
		if(mark == '5'){ *(*(mi+0)+0) = 1;*(*(mi+1)+0) = 0;}//aa
                else if(mark == '6'){*(*(mi+0)+0) = 0;*(*(mi+1)+0) = 1;}//ab
                else{
			//fprintf(stderr,"Error: aaxab phenotype mark's type must be(aa,ab).\n ");
		}
	}else {
		fprintf(stderr,"phenotype error.\n");
		exit(0);
	}
}

void get_H(double r,double *H_point,int phase){
	//phase暂时没用
	double (*H)[4] = (double (*)[4])H_point;
	
	if(phase ==  1){
		*(*(H+0)+0) = (1 - r) * (1 - r);	*(*(H+0)+1) = r * (1 - r);		*(*(H+0)+2) = r * (1 - r);		*(*(H+0)+3) = r * r;
		*(*(H+1)+0) = r * (1 - r);		*(*(H+1)+1) = (1 - r) * (1 - r);	*(*(H+1)+2) = r * r;			*(*(H+1)+3) = r * (1 - r);
		*(*(H+2)+0) = r * (1 - r);		*(*(H+2)+1) = r * r;			*(*(H+2)+2) = (1 - r) * (1 - r);	*(*(H+2)+3) = r * (1 - r);
		*(*(H+3)+0) = r * r;			*(*(H+3)+1) = r * (1 - r);		*(*(H+3)+2) = r * (1 - r);		*(*(H+3)+3) = (1 - r) * (1 - r);
	}else if(phase == 2){
		*(*(H+0)+0) = (1 - r) * r;	        *(*(H+0)+1) = (1 - r) * (1 - r);        *(*(H+0)+2) = r * r;	                *(*(H+0)+3) = (1 - r) * r;
                *(*(H+1)+0) = (1 - r) * (1 - r);        *(*(H+1)+1) = r * (1 - r); 	        *(*(H+1)+2) = (1 - r) * r;              *(*(H+1)+3) = r * r;
                *(*(H+2)+0) = r * r;  	                *(*(H+2)+1) = r * (1 - r);              *(*(H+2)+2) = (1 - r) * r;              *(*(H+2)+3) = (1 - r) * (1 - r);
                *(*(H+3)+0) = r * (1 - r);              *(*(H+3)+1) = r * r;          	        *(*(H+3)+2) = (1 - r) * (1 - r);        *(*(H+3)+3) = (1 - r) * r;
	}else if(phase == 3){
		*(*(H+0)+0) = (1 - r) * r;              *(*(H+0)+1) = r * r;                    *(*(H+0)+2) = r * r;	                *(*(H+0)+3) = (1 - r) * r;
                *(*(H+1)+0) = r * r;                    *(*(H+1)+1) = (1 - r) * r;              *(*(H+1)+2) = (1 - r) * r;              *(*(H+1)+3) = r * r;
                *(*(H+2)+0) = (1 - r) * (1 - r);        *(*(H+2)+1) = (1 -  r) * r;             *(*(H+2)+2) = (1 - r) * r;              *(*(H+2)+3) = (1 - r) * (1 - r);
                *(*(H+3)+0) = (1 - r) * r;              *(*(H+3)+1) = (1 - r) * (1 - r);        *(*(H+3)+2) = (1 - r) * (1 - r);        *(*(H+3)+3) = (1 - r) * r;
        }else if(phase == 4){
		*(*(H+0)+0) = r * r;                    *(*(H+0)+1) = r * (1 - r);              *(*(H+0)+2) = r * (1 - r);              *(*(H+0)+3) = (1 - r) * (1 - r);
                *(*(H+1)+0) = r * (1 - r);              *(*(H+1)+1) = r * r;                    *(*(H+1)+2) = (1 -  r) * (1 - r);        *(*(H+1)+3) = r * (1 - r);
                *(*(H+2)+0) = r * (1 - r);              *(*(H+2)+1) = (1 - r) * (1 - r);        *(*(H+2)+2) = r * r;                     *(*(H+2)+3) = r * (1 - r);
                *(*(H+3)+0) = (1 - r) * (1 - r);        *(*(H+3)+1) = r * (1 - r);              *(*(H+3)+2) = r * (1 - r);               *(*(H+3)+3) = r * r;
        }
}

void get_D(double *D_point,int phase){
	//同样 pahse暂时没用
	double (*D)[4] = (double (*)[4])D_point;

	if(phase == 1){
		*(*(D+0)+0) = 0; *(*(D+0)+1) = 1; *(*(D+0)+2) = 1; *(*(D+0)+3) = 2;
		*(*(D+1)+0) = 1; *(*(D+1)+1) = 0; *(*(D+1)+2) = 2; *(*(D+1)+3) = 1;
		*(*(D+2)+0) = 1; *(*(D+2)+1) = 2; *(*(D+2)+2) = 0; *(*(D+2)+3) = 1;
		*(*(D+3)+0) = 2; *(*(D+3)+1) = 1; *(*(D+3)+2) = 1; *(*(D+3)+3) = 0;
	}else if(phase == 2){
		*(*(D+0)+0) = 1; *(*(D+0)+1) = 0; *(*(D+0)+2) = 2; *(*(D+0)+3) = 1;
                *(*(D+1)+0) = 0; *(*(D+1)+1) = 1; *(*(D+1)+2) = 1; *(*(D+1)+3) = 2;
                *(*(D+2)+0) = 2; *(*(D+2)+1) = 1; *(*(D+2)+2) = 1; *(*(D+2)+3) = 0;
                *(*(D+3)+0) = 1; *(*(D+3)+1) = 2; *(*(D+3)+2) = 0; *(*(D+3)+3) = 1;
	}else if(phase == 3){
		*(*(D+0)+0) = 1; *(*(D+0)+1) = 2; *(*(D+0)+2) = 0; *(*(D+0)+3) = 1;
                *(*(D+1)+0) = 2; *(*(D+1)+1) = 1; *(*(D+1)+2) = 1; *(*(D+1)+3) = 0;
                *(*(D+2)+0) = 0; *(*(D+2)+1) = 1; *(*(D+2)+2) = 1; *(*(D+2)+3) = 2;
                *(*(D+3)+0) = 1; *(*(D+3)+1) = 0; *(*(D+3)+2) = 2; *(*(D+3)+3) = 1;
        }else if(phase == 4){
		*(*(D+0)+0) = 2; *(*(D+0)+1) = 1; *(*(D+0)+2) = 1; *(*(D+0)+3) = 0;
                *(*(D+1)+0) = 1; *(*(D+1)+1) = 2; *(*(D+1)+2) = 0; *(*(D+1)+3) = 1;
                *(*(D+2)+0) = 1; *(*(D+2)+1) = 0; *(*(D+2)+2) = 2; *(*(D+2)+3) = 1;
                *(*(D+3)+0) = 0; *(*(D+3)+1) = 1; *(*(D+3)+2) = 1; *(*(D+3)+3) = 2;
        }
}	

int get_I_row(int phenotype){
	if(phenotype == 1 || phenotype == 2){
		return 4;
	}else if(phenotype == 3){
		return 3;
	}else if(phenotype == 4 || phenotype == 5){
		return 2;
	}else{
		fprintf(stderr,"phase error.\n");
		exit(0);
	}
}

void get_I_array(double *I_point,int row,int phenotype){
	double (*I)[row] = (double (*)[row])I_point;
	if(phenotype == 1 || phenotype == 2){
		*(*(I+0)+0) = 1;*(*(I+0)+1) = 0;*(*(I+0)+2) = 0;*(*(I+0)+3) = 0;
		*(*(I+1)+0) = 0;*(*(I+1)+1) = 1;*(*(I+1)+2) = 0;*(*(I+1)+3) = 0;
		*(*(I+2)+0) = 0;*(*(I+2)+1) = 0;*(*(I+2)+2) = 1;*(*(I+2)+3) = 0;
		*(*(I+3)+0) = 0;*(*(I+3)+1) = 0;*(*(I+3)+2) = 0;*(*(I+3)+3) = 1;
	}else if(phenotype == 3){
		*(*(I+0)+0) = 1;*(*(I+0)+1) = 0;*(*(I+0)+2) = 0;
		*(*(I+1)+0) = 0;*(*(I+1)+1) = 1;*(*(I+1)+2) = 0;
		*(*(I+2)+0) = 0;*(*(I+2)+1) = 1;*(*(I+2)+2) = 0;
		*(*(I+3)+0) = 0;*(*(I+3)+1) = 0;*(*(I+3)+2) = 1;
	}else if(phenotype == 4){
		*(*(I+0)+0) = 1;*(*(I+0)+1) = 0;
		*(*(I+1)+0) = 1;*(*(I+1)+1) = 0;
		*(*(I+2)+0) = 0;*(*(I+2)+1) = 1;
		*(*(I+3)+0) = 0;*(*(I+3)+1) = 1;
	}else if(phenotype == 5){
		*(*(I+0)+0) = 1;*(*(I+0)+1) = 0;
                *(*(I+1)+0) = 0;*(*(I+1)+1) = 1;
                *(*(I+2)+0) = 1;*(*(I+2)+1) = 0;
                *(*(I+3)+0) = 0;*(*(I+3)+1) = 1;
	}else{
		fprintf(stderr,"phenotype error.\n");exit(0);
	}
}


//cp所有群体类型
void cp(int phenotype1,int phenotype2,char *mark1_array,char *mark2_array,double *result,double r0,double accuracy,int SAMP,int LOD_flag){

	//LOD_flag == 1时，只做LRT和LOD计算
	if((!LOD_flag) && phenotype1 + phenotype2 == 9){//nnxnp 和 lmxll类型无法计算重组率
		*(result + 0) = 1;
		*(result + 3) = 0;
		return;
	}
	double temp_MLE[4];
	double temp_LRT[4];
	double temp_LOD[4];
	int phase = 1;//两个mark的遗传相
	while(phase <= 4){
		double distance = 1;
		double r = r0;
		int final_flag = 1;
		int ln_flag = LOD_flag;//用于计算LOD.为了代码的重用性，不在另写函数,也用于纠正
		int ln_05_flag = 0;
		double fenzi = 1;
		double fenmu = 1;
		int cishu = 0;

		do{
			cishu++;
			double temp = 0;
			int i;
			for(i=0;i<SAMP;i++){
				char mark1 = *(mark1_array + i);
				char mark2 = *(mark2_array + i);
				if(mark1 == '8' || mark2 == '8'){continue;}
				double H_array[4][4];
				double *H = (double *)H_array;
				get_H(r,H,phase);

				double D_array[4][4];
				double *D = (double *)D_array;
				get_D(D,phase);
///////////////////////////mark1/////////////////////////////
				int b1 = get_I_row(phenotype1);
	
				double mi_j1[b1][1];
				get_mi_matrix(phenotype1,mark1,mi_j1);

				double mi_j1_T[1][b1];//Mi_j1_T(1,b1)
				double *mi_j1_point = (double *)mi_j1;
				double *mi_j1_T_point = (double *)mi_j1_T;
				matrix_transp(mi_j1_point,b1,1,mi_j1_T_point,1,b1);
		
				double I_b1[4][b1];
				double *I_b1_point = (double *)I_b1;
				get_I_array(I_b1_point,b1,phenotype1);
		
				double I_b1_T[b1][4];//I_b1_T (b1,4)
				double *I_b1_T_point = (double *)I_b1_T;
				matrix_transp(I_b1_point,4,b1,I_b1_T_point,b1,4);
////////////////////////////mark2/////////////////////////////
				int b2 = get_I_row(phenotype2);
	
				double mi_j2[b2][1];
				double *mi_j2_point = (double *)mi_j2;
				get_mi_matrix(phenotype2,mark2,mi_j2);

				double I_b2[4][b2];
                		double *I_b2_point = (double *)I_b2;
		                get_I_array(I_b2_point,b2,phenotype2);
////////////////////////////////////////////////////////////
				double HD[4][4];
				double *HD_point = (double *)HD;
				matrix_elementwise(H,4,4,D,4,4,HD_point,4,4);

				double matrix1[b1][4];
				double *matrix1_point = (double *)matrix1;
				matrix_multi(I_b1_T_point,b1,4,HD_point,4,4,matrix1_point,b1,4);

				double matrix2[b1][b2];
				double *matrix2_point = (double *)matrix2;
				matrix_multi(matrix1_point,b1,4,I_b2_point,4,b2,matrix2_point,b1,b2);

				double matrix3[1][b2];
				double *matrix3_point = (double *)matrix3;
				matrix_multi(mi_j1_T_point,1,b1,matrix2_point,b1,b2,matrix3_point,1,b2);

				double matrix4[1][1];
				double *matrix4_point = (double *)matrix4;
				matrix_multi(matrix3_point,1,b2,mi_j2_point,b2,1,matrix4_point,1,1);
//--------------------------------------------------------------------------------------以上是MLE分子，以下是MLE分母
				double matrix5[b1][4];
                		double *matrix5_point = (double *)matrix5;
		                matrix_multi(I_b1_T_point,b1,4,H,4,4,matrix5_point,b1,4);

                		double matrix6[b1][b2];
		                double *matrix6_point = (double *)matrix6;
                		matrix_multi(matrix5_point,b1,4,I_b2_point,4,b2,matrix6_point,b1,b2);

		                double matrix7[1][b2];
                		double *matrix7_point = (double *)matrix7;
		                matrix_multi(mi_j1_T_point,1,b1,matrix6_point,b1,b2,matrix7_point,1,b2);

                		double matrix8[1][1];
		                double *matrix8_point = (double *)matrix8;
                		matrix_multi(matrix7_point,1,b2,mi_j2_point,b2,1,matrix8_point,1,1);
//---------------------------------------------------------------------------------------
				temp += matrix4[0][0] / matrix8[0][0];
	
				if(ln_flag && (!ln_05_flag)){fenzi *= matrix8[0][0];}//计算Ln(r)
				if(ln_05_flag){fenmu *= matrix8[0][0];}//计算Ln(r)
			}
			double new_r = temp / SAMP / 2;
			distance = fabs(new_r - r);
			//if(!ln_flag){printf("new_r=%lf\tdistance=%lf\n",new_r,distance);}
			r = new_r;

			if(ln_05_flag){final_flag = 0;}//LOD分母计算完毕，全部执行完毕，退出                                              /^\ 激活顺序向上
			if(ln_flag){ln_05_flag = 1;r = 0.5;}//LOD分子计算完毕，激活LOD分母计算,将R指改成0.5	                           |
			if((!ln_flag) && ((distance <= accuracy) || (cishu > 100))){ln_flag = 1;temp_MLE[phase - 1] = r;}  //重组率计算完毕，激活LOD分子计算
	
		}while(final_flag);

		temp_LRT[phase - 1] = (-2) * log(fenzi / fenmu);
		temp_LOD[phase - 1] =  log10(fenzi / fenmu);
		//printf("**********%lf\t%lf\n",temp_LRT[phase - 1],temp_LOD[phase - 1]);
		phase++;
	}
//--------------------------------选出最好的重组率,先选出一个小于0.5. 然后根据LOD越大越准确的原则
	int best_index = 0;
	while(best_index < 4){
		if(temp_MLE[best_index] < 0.5){
			break;
		}else{
			best_index++;
		}
	}

	int temp_index = 0;
	while( temp_index < 4){
		if(temp_MLE[temp_index] < 0.5 && temp_LOD[temp_index] >= temp_LOD[best_index]){
			best_index = temp_index;
		}
		temp_index++;
	}

	*(result + 0) = temp_MLE[best_index];
	*(result + 1) = 0;
	*(result + 2) = temp_LRT[best_index];
	*(result + 3) = temp_LOD[best_index];
}
