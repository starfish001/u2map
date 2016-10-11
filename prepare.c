#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calculate.h"
#include "correct.h"
#include "prepare.h"
int step1_prepare(char *inputFileName,char *outputFileName,int grouptype,float accuracy){

	//参数设置
	int MARK_NUM = 0;
	int SAMP_NUM = 0;
	//printf("%f\n",accuracy);
	if(accuracy == 0){accuracy = 0.001;}
	double ACCURACY = accuracy;
	double R = 0.25;
	
	//判断程序是否正确调用
	if(!(inputFileName && outputFileName && grouptype))errorPrint();

	if(grouptype < 1 || grouptype > 21)errorPrint();

	FILE *infile;//输入文件
	if(!(infile = fopen(inputFileName,"r"))){
		fprintf(stderr,"Cannot open file data.txt, please check.\n");
                exit(0);
	}

	FILE *tempfile;//临时输出文件,包含了MLE,VAR,LRT,LOD所有信息

	int len = strlen(outputFileName);
	char TempName[len + 12];
	strcpy(TempName,outputFileName);
	strcat(TempName,"_temp.list");
        if(!(tempfile = fopen(TempName,"w"))){
                fprintf(stderr,"Cannot open temp file,output directory not exist. please check.\n");
                exit(0);
        }
	fputs("MLE\tVAR\tLRT\tLOD\n",tempfile);

	//读取mark数和个体数
	fscanf(infile,"MARK_NUM=%d\n",&MARK_NUM);
	fscanf(infile,"INDVIDUAL_NUM=%d\n",&SAMP_NUM);
	if(SAMP_NUM == 0 && MARK_NUM == 0){fprintf(stderr,"Read file error,input file must begin with MARK_NUM and INDVIDUAL_NUM\n");exit(EXIT_FAILURE);}

	/*去除第一行输入，（20151113 陈海新）
	跳过表格第一行,CP群体也一样
	char header;
	header = fgetc(infile); 
	int idiv = 0;
	while(infile){
		header = fgetc(infile);
		if(header == '\t')idiv++;//注意所有样品的个数必须和INDVIDUAL_NUM保持一致
		if(header == '\n')break;
	}

	if(idiv > SAMP_NUM){
		//samples number of inputfile more then SAMP_NUM
		fprintf(stderr,"Warn: The samples number is more then %d ,but just analyse %d samples .\n",SAMP_NUM,SAMP_NUM);
		SAMP_NUM =  idiv;
	}else if(idiv < SAMP_NUM){
		//samples number of inputfile less then SAMP_NUM
		fprintf(stderr,"Error: The samples number is less then %d , please check. \n",SAMP_NUM);
                exit(0);
	}
	*/
	int i,j;
	long offset = 0;//offset for fseek
	char markName[MARK_NUM][30];//所哟mark的名称 

	char mark1[SAMP_NUM];//mark1 of all sample
	char mark2[SAMP_NUM];//mark2 of all sample
	char *mark_p1 = mark1;//mark1 point of all sample
	char *mark_p2 = mark2;//mark2 point of all sample

	double result[4];//calculate result of two mark
	double *result_p = result;//result point
	float percentage1[3];
	float percentage2[3];//N的可能值概率，percentage2[0] 是P1，percentage2[1] 是P2，percentage2[2] 是F1
	float* percentage1_p = percentage1;
	float* percentage2_p = percentage2;
	//计算各值，并输出到临时文件
	offset = ftell(infile);
	for(i = 0;i < MARK_NUM;i++){
		char temp_name1[30] = {"\0"};//mark1 名称

		int maketype1,maketype2;
		fseek(infile,offset,0);
		maketype1 = readline(infile,mark_p1,temp_name1,SAMP_NUM,percentage1_p,grouptype);
		strcpy(markName[i],temp_name1);
		//printf("%d,%s\n",i,markName[i]);

		offset = ftell(infile);
		for(j = i + 1;j < MARK_NUM;j++){
			char temp_name2[30] = {"\0"};//mark2 名称
			//printf("%d\t",j);
			result[0] = result[1] = result[2] = result[3] = 0;
			maketype2 = readline(infile,mark_p2,temp_name2,SAMP_NUM,percentage2_p,grouptype);

			if(grouptype <= 20){//非CP群体
				float  eachN[9] = {0};//n1,n2,...n9 (AABB,aabb,AaBb,AaBB...)
				//printf("%d,%d\t",maketype1,maketype2);
				statMark(temp_name1,temp_name2,mark_p1,mark_p2,eachN,SAMP_NUM,percentage1,percentage2,grouptype);
				calculate(ACCURACY,R,eachN,result_p,grouptype,maketype1,maketype2);
				if(result[0] == 0)
				{
					result[3] = 50;//重组率为零时,计算的LOD为nan值
				}else if(result[0] < 0 || result[0] > 0.5)//无法准确计算重组率，则折中处理
				{
					result[0] = 0.25;
					result[3] = 1;
				}else
				{
					result[3] = calLOD(result[2],result[0]);//LOD,CP群体已计算好LOD
				}
				//printf("finish!\n");

				fprintf(tempfile,"%s\t%s\t%lf\t%lf\t%lf\t%lf\n",temp_name1,temp_name2,result[0],result[1],result[2],result[3]);
			}else{//CP群体
				cp(maketype1,maketype2,mark_p1,mark_p2,result_p,R,ACCURACY,SAMP_NUM,0);
				//printf("finish.\n");
				fprintf(tempfile,"%s\t%s\t%lf\t%lf\t%lf\t%lf\n",temp_name1,temp_name2,result[0],result[1],result[2],result[3]);
				//fprintf(tempfile,"%lf\t%lf\t%lf\t%lf\n",result[0],result[1],result[2],result[3]);
				//printf("%s\t%s\t%lf\t%lf\t%lf\t%lf\n",temp_name1,temp_name2,result[0],result[1],result[2],result[3]);
			}
		}
	}
	fclose(infile);
	fclose(tempfile);


	//临时文件转化成四个MLE，VAR，LRT，LOD
	FILE *MLE,*VAR,*LRT,*LOD;
	char MLENAME[len + 15];
	char VARNAME[len + 15];
	char LRTNAME[len + 15];
	char LODNAME[len + 15];
	strcpy(MLENAME,outputFileName);
	strcpy(VARNAME,outputFileName);
	strcpy(LRTNAME,outputFileName);
	strcpy(LODNAME,outputFileName);
	strcat(MLENAME,".MLE.table.xls");
	strcat(VARNAME,".VAR.table.xls");
	strcat(LRTNAME,".LRT.table.xls");
	strcat(LODNAME,".LOD.table.xls");
	
	if(!(tempfile = fopen(TempName,"r"))){fprintf(stderr,"Cannot open file data.txt, please check.\n");exit(0);}
	if(!(MLE = fopen(MLENAME,"w"))){fprintf(stderr,"Cannot open MLE.table.xls, please check.\n");exit(0);}
	if(!(VAR = fopen(VARNAME,"w"))){fprintf(stderr,"Cannot open VAR.table.xls, please check.\n");exit(0);}
	if(!(LRT = fopen(LRTNAME,"w"))){fprintf(stderr,"Cannot open LRT.table.xls, please check.\n");exit(0);}
	if(!(LOD = fopen(LODNAME,"w"))){fprintf(stderr,"Cannot open LOD.table.xls, please check.\n");exit(0);}

	fprintf(MLE,"markNum=%d\n",MARK_NUM);
	fprintf(VAR,"markNum=%d\n",MARK_NUM);
	fprintf(LRT,"markNum=%d\n",MARK_NUM);
	fprintf(LOD,"markNum=%d\n",MARK_NUM);
	
	char title[30];
	float data;
	fscanf(tempfile,"%s",title);//跳过MLE
	fscanf(tempfile,"%s",title);//跳过VAR
	fscanf(tempfile,"%s",title);//跳过LRT
	fscanf(tempfile,"%s",title);//跳过LOD
	i = 1;j = 1;
	fprintf(MLE,"%s",markName[0]);
	fprintf(VAR,"%s",markName[0]);
	fprintf(LRT,"%s",markName[0]);
	fprintf(LOD,"%s",markName[0]);
	while(!feof(tempfile))
	{
		i++;
		fscanf(tempfile,"%s",title);//跳过mark1
		fscanf(tempfile,"%s",title);//跳过mark2
		fscanf(tempfile,"%f",&data);
		fprintf(MLE,"\t%f",data);
		fscanf(tempfile,"%f",&data);
		fprintf(VAR,"\t%f",data);
		fscanf(tempfile,"%f",&data);
		fprintf(LRT,"\t%f",data);
		fscanf(tempfile,"%f\n",&data);
		fprintf(LOD,"\t%f",data);
		if(i == MARK_NUM)
		{
			fprintf(MLE,"\n%s",markName[j]);
			fprintf(VAR,"\n%s",markName[j]);
			fprintf(LRT,"\n%s",markName[j]);
			fprintf(LOD,"\n%s",markName[j]);
			j++;
			i = j;
		}
	}

	fclose(tempfile);
	fclose(MLE);
	fclose(VAR);
	fclose(LRT);
	fclose(LOD);
	remove(TempName);//删除临时文件

	//CP群体对结果进行纠正
	if(grouptype == 21){
		fprintf(stdout,"CP group, correcting result...\n");
		correct_CP_MLE(inputFileName,outputFileName,accuracy);
		fprintf(stdout,"correcting finish!\n");
	}
	
	return 0;

}

//calculate MLE VAR LRT
void calculate(double accuracy,double r,float eachN[],double *result,int grouptype,int marktype1,int marktype2){	
	switch(grouptype){
		case 1 : f1dh(eachN,result);break;
		case 2 : f1ril(eachN,result);break;
		case 3 : p1bc1ril(eachN,result,r,accuracy);break;
		case 4 : p2bc1ril(eachN,result,r,accuracy);break;
		case 5 : p1bc2ril(eachN,result,r,accuracy);break;
		case 6 : p2bc2ril(eachN,result,r,accuracy);break;
		case 7 : p1bc1dh(eachN,result,r,accuracy);break;
		case 8 : p2bc1dh(eachN,result,r,accuracy);break;
		case 9 : p1bc2dh(eachN,result,r,accuracy);break;
		case 10: p2bc2dh(eachN,result,r,accuracy);break;
		case 11: p1bc1f1(eachN,result,r,accuracy);break;
		case 12: p2bc1f1(eachN,result,r,accuracy);break;
		case 13: f2(eachN,result,r,accuracy);break;
		case 14: f3g(eachN,result,r,accuracy);break;
		case 15: p1bc2f1(eachN,result,r,accuracy);break;
		case 16: p2bc2f1(eachN,result,r,accuracy);break;
		case 17: p1bc1f2(eachN,result,r,accuracy);break;
		case 18: p2bc1f2(eachN,result,r,accuracy);break;
		case 19: p1bc2f2(eachN,result,r,accuracy);break;
		case 20: p2bc2f2(eachN,result,r,accuracy);break;
		default : fprintf(stderr,"please check group number.\n");
	}
}

//read one line of sample_mark file
int readline(FILE * fp,char *mark,char *mark_name,int SAMP,float *percentage,int grouptype){
	int marktype = 0;//0代表双纯，1代表lm*ll，2代表nn*np，3代表双杂合
	char type[10];
	if(grouptype <= 20){
	        fscanf(fp,"%s",mark_name);
		//fscanf(fp,"%f:%f:%f",percentage+0,percentage+1,percentage+2);
                	int idiv = 0;
	                char temp;
        	        do{
                	        temp = fgetc(fp);
                        	if(temp != '\t'){
                                	*(mark + idiv) = temp;
	                                idiv++;
        	                }
                	        if(idiv >= SAMP){//if samples number is more then SAMP_NUM,skip extra samples.
                        	        do{
                                	        temp = fgetc(fp);
	                                }while(temp != '\n');
        	                }
                	}while(temp != '\n');
			if(idiv < SAMP){fprintf(stderr,"ERROR: %s's gene type is less than %d.\n",mark_name,SAMP);exit(EXIT_FAILURE);}
	}else{//CP群体,20150811添加。陈海新
		fscanf(fp,"%s",mark_name);
		fscanf(fp,"%s",type);
		//printf("%s\t%s...",mark_name,type);
		if(strcmp(type,"abxcd") == 0){marktype = 1;}
		else if(strcmp(type,"abxac") == 0){marktype = 2;}
		else if(strcmp(type,"abxab") == 0){marktype = 3;}
		else if(strcmp(type,"abxaa") == 0){marktype = 4;}
		else if(strcmp(type,"aaxab") == 0){marktype = 5;}
		else{fprintf(stderr,"ERROR: %s's gene type error,can not be %s.\n",mark_name,type);exit(EXIT_FAILURE);}
                //fscanf(fp,"%f:%f:%f",percentage+0,percentage+1,percentage+2);
                        int idiv = 0;
                        char temp_type[4];
			//segment  = fgetc(fp);
                        do{
				fscanf(fp,"%s",temp_type);
				if(strcmp(temp_type,"ac") == 0 || strcmp(temp_type,"ca") == 0){*(mark + idiv) = '1';}
				else if(strcmp(temp_type,"ad") == 0 || strcmp(temp_type,"da") == 0){*(mark + idiv) = '2';}
				else if(strcmp(temp_type,"bc") == 0 || strcmp(temp_type,"cb") == 0){*(mark + idiv) = '3';}
				else if(strcmp(temp_type,"bd") == 0 || strcmp(temp_type,"db") == 0){*(mark + idiv) = '4';}
				else if(strcmp(temp_type,"ab") == 0 || strcmp(temp_type,"ba") == 0){*(mark + idiv) = '6';}
				else if(strcmp(temp_type,"aa") == 0){*(mark + idiv) = '5';}
				else if(strcmp(temp_type,"bb") == 0){*(mark + idiv) = '7';}
				else{
					*(mark + idiv) = '8';
					//printf("hello %s\n",temp_type);
				}

                                /*if(idiv >= SAMP){//if samples number is more then SAMP_NUM,skip extra samples.
                                        do{
                                        	segment = fgetc(fp);
                                        }while(segment != '\n');
                                }*/

                                idiv++;
                        }while(idiv < SAMP);
			if(idiv < SAMP){fprintf(stderr,"ERROR: %s's gene type is less than %d.\n",mark_name,SAMP);exit(EXIT_FAILURE);}
	}
	return marktype;
}

//calculate LOD
double calLOD(double lrt,double mle){
	if(mle < 0 ){//重组率小于0
		return 0;
	}else{
	        return lrt/(2*log(10));
	}
}


//statistics N1,N2,N3,N4 of two marks
void statMark(char *make_name1,char *make_name2,char *mark1,char *mark2,float *eachN,int SAMP,float percentage1[],float percentage2[],int grouptype){
        int i;
        char character1,character2;

	if(grouptype < 21){
	        for(i = 0; i < SAMP; i++){
        	        character1 = *(mark1 + i);
	                character2 = *(mark2 + i);
	       	        if(character1 == 'X'){
        	                if(character2 == 'X'){
                	                *(eachN + 0) = *(eachN + 0) + 1;
                        	}else if(character2 == 'Z'){
	                                *(eachN + 1) = *(eachN + 1) + 1;
        	                }else if(character2 == 'Y'){
					*(eachN + 2) = *(eachN + 2) + 1;
				}else if(character2 == 'N'){
				//*(eachN + 0) = *(eachN + 0) + percentage2[0];
				//*(eachN + 1) = *(eachN + 1) + percentage2[1];
				//*(eachN + 2) = *(eachN + 2) + percentage2[2];
					continue;
				}else{
					fprintf(stderr,"ERROR,%s gene type wrong ,can not be %c.\n",make_name2,character2);exit(EXIT_FAILURE);
				}
	                }else if(character1 == 'Z'){
        	                if(character2 == 'X'){
                	                *(eachN + 3) = *(eachN + 3) + 1;
	                        }else if(character2 == 'Z'){
        	                        *(eachN + 4) = *(eachN + 4) + 1;
                	        }else if(character2 == 'Y'){
					*(eachN + 5) = *(eachN + 5) + 1;
				}else if(character2 == 'N'){
				//*(eachN + 3) = *(eachN + 3) + percentage2[0];
				//*(eachN + 4) = *(eachN + 4) + percentage2[1];
				//*(eachN + 5) = *(eachN + 5) + percentage2[2];
					continue;
				}else{
                                        fprintf(stderr,"ERROR,%s gene type wrong ,can not be %c.\n",make_name2,character2);exit(EXIT_FAILURE);
                                }
	               	}else if(character1 == 'Y'){
				if(character2 == 'X'){
					*(eachN + 6) = *(eachN + 6) + 1;
				}else if(character2 == 'Z'){
					*(eachN + 7) = *(eachN + 7) + 1;
				}else if(character2 == 'Y'){
					*(eachN + 8) = *(eachN + 8) + 1;
				}else if(character2 == 'N'){
				//*(eachN + 6) = *(eachN + 6) + percentage2[0];
				//*(eachN + 7) = *(eachN + 7) + percentage2[1];
				//*(eachN + 8) = *(eachN + 8) + percentage2[2];
					continue;
				}
			}else if(character1 == 'N'){
					continue;
/*
			if(character2 == 'X'){
				*(eachN + 0) = *(eachN + 0) + percentage1[0];
				*(eachN + 3) = *(eachN + 3) + percentage1[1];
				*(eachN + 6) = *(eachN + 6) + percentage1[2];
			}else if(character2 == 'Z'){
				*(eachN + 1) = *(eachN + 1) + percentage1[0];
				*(eachN + 4) = *(eachN + 4) + percentage1[1];
				*(eachN + 7) = *(eachN + 7) + percentage1[2];
			}else if(character2 == 'Y'){
				*(eachN + 2) = *(eachN + 2) + percentage1[0];
				*(eachN + 5) = *(eachN + 5) + percentage1[1];
				*(eachN + 8) = *(eachN + 8) + percentage1[2];
			}else if(character2 == 'N'){
				*(eachN + 0) = *(eachN + 0) + percentage1[0] * percentage2[0];
				*(eachN + 1) = *(eachN + 1) + percentage1[0] * percentage2[1];
				*(eachN + 2) = *(eachN + 2) + percentage1[0] * percentage2[2];
				*(eachN + 3) = *(eachN + 3) + percentage1[1] * percentage2[0];
                                *(eachN + 4) = *(eachN + 4) + percentage1[1] * percentage2[1];
                                *(eachN + 5) = *(eachN + 5) + percentage1[1] * percentage2[2];
                                *(eachN + 6) = *(eachN + 6) + percentage1[2] * percentage2[0];
                                *(eachN + 7) = *(eachN + 7) + percentage1[2] * percentage2[1];
                                *(eachN + 8) = *(eachN + 8) + percentage1[2] * percentage2[2];

			}*/
			}else{
                                fprintf(stderr,"ERROR,%s gene type wrong ,can not be %c.\n",make_name1,character1);exit(EXIT_FAILURE);
                        }
		}
        }else{//CP 20150817 陈海新
		for(i = 0; i < SAMP; i++){
                        character1 = *(mark1 + i);
                        character2 = *(mark2 + i);
                        if(character1 == 'A'){
                                if(character2 == 'A'){
                                        *(eachN + 0) = *(eachN + 0) + 1;//AABB
                                }else if(character2 == 'Z'){
                                        *(eachN + 1) = *(eachN + 1) + 1;//AABb
                                }else if(character2 == 'a'){
                                        *(eachN + 2) = *(eachN + 2) + 1;//AAbb
                                }else if(character2 == 'N'){
                                        continue;
                                }else{
                                        fprintf(stderr,"ERROR,%s gene type wrong ,can not be %c.\n",make_name2,character2);exit(EXIT_FAILURE);
                                }
                        }else if(character1 == 'Z'){
                                if(character2 == 'A'){
                                        *(eachN + 3) = *(eachN + 3) + 1;//AaBB
                                }else if(character2 == 'Z'){
                                        *(eachN + 4) = *(eachN + 4) + 1;//AaBb
                                }else if(character2 == 'a'){
                                        *(eachN + 5) = *(eachN + 5) + 1;//Aabb
                                }else if(character2 == 'N'){
                                        continue;
                                }else{
                                        fprintf(stderr,"ERROR,%s gene type wrong ,can not be %c.\n",make_name2,character2);exit(EXIT_FAILURE);
                                }
                        }else if(character1 == 'a'){
                                if(character2 == 'A'){
                                        *(eachN + 6) = *(eachN + 6) + 1;//aaBB
                                }else if(character2 == 'Z'){
                                        *(eachN + 7) = *(eachN + 7) + 1;//aaBb
                                }else if(character2 == 'a'){
                                        *(eachN + 8) = *(eachN + 8) + 1;//aabb
                                }else if(character2 == 'N'){
                                        continue;
                                }else{
                                        fprintf(stderr,"ERROR,%s gene type wrong ,can not be %c.\n",make_name2,character2);exit(EXIT_FAILURE);
                                }
                        }else if(character1 == 'N'){
                                continue;
			}else{
                                fprintf(stderr,"ERROR,%s gene type wrong ,can not be %c.\n",make_name1,character1);exit(EXIT_FAILURE);
                        }
		}
	}
	//int x;
	//for(x = 0;x<9;x++){printf("%f\t",*(eachN + x));}
	//printf("\n");
}

//软件输入形式
void errorPrint()
{
    fprintf(stderr,"Usage:u2map prepare -f1 mark.input -o . -g 1 -p prefix -a 0.001\n");
    fprintf(stderr,"\n\t-f1\t[str] Input file,begin with markNum.\n");
    fprintf(stderr,"\t-o\t[str] Output directory,must exist. [default .]\n");
    fprintf(stderr,"\t-p\t[str] Output file prefix. [default prefix]\n");
    fprintf(stderr,"\t-a\t[float] Accuracy of recombination frequency. [default 0.001]\n");
    fprintf(stderr,"\t-g\t[int] Group type, from 1 to 20.\n");
    fprintf(stderr,"\t\t 01:F1DH      02:F1RIL     03:P1BC1RIL \n");
    fprintf(stderr,"\t\t 04:P2BC1RIL  05:P1BC2RIL  06:P2BC2RIL \n");
    fprintf(stderr,"\t\t 07:P1BC1DH   08:P2BC1DH   09:P1BC2DH  \n");
    fprintf(stderr,"\t\t 10:P2BC2DH   11:P1BC1F1   12:P2BC1F1  \n");
    fprintf(stderr,"\t\t 13:F2        14:F3        15:P1BC2F1  \n");
    fprintf(stderr,"\t\t 16:P2BC2F1   17:P1BC1F2   18:P2BC1F2  \n");
    fprintf(stderr,"\t\t 19:P1BC2F2   20:P2BC2F2   21.CP   \n\n");

    exit(EXIT_FAILURE);

}

