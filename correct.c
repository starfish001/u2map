//用于mle纠正，只针对CP群体
//CP群体nnxnp和lmxll类型的mark无法计算重组率，计算时赋值为1;其他群体重组率大于0.5的marks也可以调整.该函数是通过其他两个mark的重组率，估算这两个mark的重组率
//为了后续的维护，同时也将方差和LOD进行纠正,因此需要读入mark的信息

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "calculate.h"
#include "correct.h"

int readline2(FILE * fp,char *mark,char *mark_name,int SAMP,int grouptype){
        int marktype = 0;//0代表双纯，1代表lm*ll，2代表nn*np，3代表双杂合
        char type[10];

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
                idiv++;
        }while(idiv < SAMP);
        
	if(idiv < SAMP){fprintf(stderr,"ERROR: %s's gene type is less than %d.\n",mark_name,SAMP);exit(EXIT_FAILURE);}
        return marktype;
}

//通过读文件获取矩阵表
void readTable2(char file[],double *r){
    FILE *fp;
    fp = fopen(file,"r");
    int markNum = 0;//MARK的数量

    fscanf(fp,"markNum=%d",&markNum);

    double (*p)[markNum] = (double (*)[markNum])r;

    int i,j;
    char name[30];
    for(i = 0;i < markNum;i++)
    {
        fscanf(fp,"%s",name);//读取MARK name
        for(j = i;j < markNum;j++)
        {
		if(i == j){
			*(*(p+i)+j) = 0;
		}else{
	                double temp;
        	        fscanf(fp,"%lf",&temp);
                	*(*(p+i)+j) = temp;
	                *(*(p+j)+i) = temp;
		}
         }
    }
    fclose(fp);
}

//重组率转为遗传距离,单位为cM
double mle2kosambi(double r){
        return 25 * log((1 + 2 * r) / (1 - 2 * r));
}

//遗传距离转化为重组率
double kosambi2mle(double k){
        double temp = exp(k/25);
        return 0.5 * ((temp - 1) / (temp + 1));
}

//判断markX 在markA和markB之间的距离,返回值为1，代表在A段；返回值为2，代表在B段，返回值为3，代表在AB的中间
int locate(double distAX,double distBX,double distAB){
        if(distAX > distAB){
                if(distAX > distBX){
                        return 2;//在B端
                }else{
                        return 1;//在A端
                }
        }else{
                if(distBX > distAB){
                        return 1;
                }else{
                        return 3;//在中间
                }
        }
}

//根据5个重组率，估计另一个重组率
double new_MLE(double distAX,double distAY,double distBX,double distBY,double distAB){
        int locateX = locate(distAX,distBX,distAB);
        int locateY = locate(distBX,distBY,distAB);

        double kosambiAX = mle2kosambi(distAX);
        double kosambiAY = mle2kosambi(distAY);
        double kosambiAB = mle2kosambi(distAB);
        double kosambiBX = mle2kosambi(distBX);
        double kosambiBY = mle2kosambi(distBY);
        double kosambiXY;

        if(locateX == 1){
                if(locateY == 1){//XY都在A端
                        kosambiXY = kosambiAX > kosambiAY ? kosambiAX - kosambiAY : kosambiAY - kosambiAX;
                }else if(locateY == 2){//X在A端，Y在B端
                        kosambiXY = kosambiAX + kosambiAB + kosambiBY;
                }else{//X在A端，Y在中间
                        kosambiXY = kosambiAX + kosambiAY;
                }
        }else if(locateX == 2){
                if(locateY == 1){
                        kosambiXY = kosambiAY + kosambiAB + kosambiBX;//X在B端，Y在A端
                }else if(locateY == 2){
                        kosambiXY = kosambiBY > kosambiBX ? kosambiBY - kosambiBX : kosambiBX - kosambiBY;//XY都在B端
                }else{
                        kosambiXY = kosambiBX + kosambiBY;//X在B端，Y在中间
                }
        }else{
                if(locateY == 1){
                        kosambiXY = kosambiAX + kosambiAY;//X在中间，Y在A端
                }else if(locateY == 2){
                        kosambiXY = kosambiBX + kosambiBY;//X在中间，Y在B端
                }else{
                        kosambiXY = kosambiAX > kosambiAY ? kosambiAX + kosambiBY - kosambiAB : kosambiAY + kosambiBX - kosambiAB;//XY在中间
                }
        }

        return kosambi2mle(kosambiXY);
}

//纠正主函数
void correct_CP_MLE(char *markFile,char *inputPrefix,double accuracy){

	int MARK_NUM = 0;
	int SAMP_NUM = 0;
//-------------------------------------------------------------读取markfile
	FILE *infile;//输入文件
        if(!(infile = fopen(markFile,"r"))){
                fprintf(stderr,"Cannot open file data.txt, please check.\n");
                exit(0);
        }
	fscanf(infile,"MARK_NUM=%d\n",&MARK_NUM);
        fscanf(infile,"INDVIDUAL_NUM=%d\n",&SAMP_NUM);
        if(SAMP_NUM == 0 && MARK_NUM == 0){fprintf(stderr,"Read file error,input file must begin with MARK_NUM and INDVIDUAL_NUM\n");exit(EXIT_FAILURE);}

        /*//跳过表格第一行,CP群体也一样
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
	char markName[MARK_NUM][30];//所有mark的名称
	char allmark[MARK_NUM][SAMP_NUM];
	char makeType[MARK_NUM];
	
	int n;
	for(n = 0;n < MARK_NUM;n++){
		makeType[n] = readline2(infile,allmark[n],markName[n],SAMP_NUM,21);
	}
	fclose(infile);
//---------------------------------------------------------------读取MLE,LRT,LOD
	int len = strlen(inputPrefix);
        char MLENAME[len + 15];
        char LRTNAME[len + 15];
        char LODNAME[len + 15];
        strcpy(MLENAME,inputPrefix);
        strcpy(LRTNAME,inputPrefix);
        strcpy(LODNAME,inputPrefix);
        strcat(MLENAME,".MLE.table.xls");
        strcat(LRTNAME,".LRT.table.xls");
        strcat(LODNAME,".LOD.table.xls");

        //double mleArray[MARK_NUM][MARK_NUM];
        //double lrtArray[MARK_NUM][MARK_NUM];
        //double lodArray[MARK_NUM][MARK_NUM];
	double (*mleArray)[MARK_NUM] = malloc(sizeof(double)*MARK_NUM*MARK_NUM);
        double (*lrtArray)[MARK_NUM] = malloc(sizeof(double)*MARK_NUM*MARK_NUM);
        double (*lodArray)[MARK_NUM] = malloc(sizeof(double)*MARK_NUM*MARK_NUM);

        double *mle = (double *)mleArray;
        double *lrt = (double *)lrtArray;
        double *lod = (double *)lodArray;

        readTable2(MLENAME,mle);
        readTable2(LRTNAME,lrt);
        readTable2(LODNAME,lod);

//-------------------------------------------------------------纠正数据，并输出到对应的文件
        FILE *newMLE,*newLRT,*newLOD;
        char newMLENAME[len + 20];
        char newLRTNAME[len + 20];
        char newLODNAME[len + 20];
        strcpy(newMLENAME,inputPrefix);
        strcpy(newLRTNAME,inputPrefix);
        strcpy(newLODNAME,inputPrefix);
        strcat(newMLENAME,".MLE.new.table.xls");
        strcat(newLRTNAME,".LRT.new.table.xls");
        strcat(newLODNAME,".LOD.new.table.xls");

        if(!(newMLE = fopen(newMLENAME,"w"))){fprintf(stderr,"Cannot open MLE.new.table.xls, please check.\n");exit(0);}
        if(!(newLRT = fopen(newLRTNAME,"w"))){fprintf(stderr,"Cannot open LRT.new.table.xls, please check.\n");exit(0);}
        if(!(newLOD = fopen(newLODNAME,"w"))){fprintf(stderr,"Cannot open LOD.new.table.xls, please check.\n");exit(0);}
//-------------开始纠正
	int x,y;
	
	fprintf(newMLE,"markNum=%d\n",MARK_NUM);
	fprintf(newLRT,"markNum=%d\n",MARK_NUM);
	fprintf(newLOD,"markNum=%d\n",MARK_NUM);
	
	for(x=0;x<MARK_NUM;x++){
		fprintf(newMLE,"%s",markName[x]);
		fprintf(newLRT,"%s",markName[x]);
		fprintf(newLOD,"%s",markName[x]);

                for(y=x;y<MARK_NUM;y++){
                        if(x==y){continue;}//X==Y,不需要输出
                        
			double temp_MLE = mleArray[x][y];
			double temp_LRT = lrtArray[x][y];
			double temp_LOD = lodArray[x][y];

                        if(temp_MLE >= 0.5){//需要调整重组率
                                int i = -1;//标记A,靠近mark X
                                int j = -1;//标记B,靠近mark Y
                                double minA = 0.5;
                                double minB = 0.5;
//先找到A，再找B，需要两次循环
                                int a = 0;
                                while(a < MARK_NUM){
                                        if(a == x || a == y){a++;continue;}
                                        double temp1 =  mleArray[a][x];
                                        double temp2 =  mleArray[a][y];
                                        if(temp1 >= 0.5 || temp2 >= 0.5){a++;continue;}//作为标杆的mark不能大于0.5
                                        if(temp1 < minA){minA = temp1; i = a;}
                                        a++;
                                }
                                a = 0;
                                while(a < MARK_NUM){
                                        if(a == x || a == y){a++;continue;}
                                        double temp1 =  mleArray[a][x];
                                        double temp2 =  mleArray[a][y];
                                        if(temp1 >= 0.5 || temp2 >= 0.5){a++;continue;}//作为标杆的mark不能大于0.5
                                        if(temp2 < minB && i != a){minB = temp2;j = a;}//i!=a;AB不能相同
                                        a++;
                                }
//找到AB的位置，开始估计
				//该位点的重组率无法估计，重组率设为0.5，LOD和LRT不变
                                if(i == -1 || j == -1){
					temp_MLE = 0.5;
					temp_LOD = 0;
					fprintf(newMLE,"\t%f",temp_MLE);
					fprintf(newLRT,"\t%f",temp_LRT);
					fprintf(newLOD,"\t%f",temp_LOD);
					continue;
				}

				//其他的可以估计
                                double distAB = mleArray[i][j];
                                double distAX = mleArray[i][x];
                                double distAY = mleArray[i][y];
                                double distBX = mleArray[j][x];
                                double distBY = mleArray[j][y];
                                double newmle;
				newmle = new_MLE(distAX,distAY,distBX,distBY,distAB);
                                if(newmle >= 0 && newmle <= 0.5){
					double result[4];
					cp(makeType[x],makeType[y],allmark[x],allmark[y],result,newmle,accuracy,SAMP_NUM,1);
					
					temp_MLE = newmle;
					temp_LRT = result[2];
					temp_LOD = result[3];
                                        fprintf(newMLE,"\t%f",temp_MLE);
                                        fprintf(newLRT,"\t%f",temp_LRT);
                                        fprintf(newLOD,"\t%f",temp_LOD);
                                }else{
                                        temp_MLE = 0.5;
					temp_LOD = 0;
					fprintf(newMLE,"\t%f",temp_MLE);
                                        fprintf(newLRT,"\t%f",temp_LRT);
                                        fprintf(newLOD,"\t%f",temp_LOD);
                                }
			 }else{
                                fprintf(newMLE,"\t%f",temp_MLE);
				fprintf(newLRT,"\t%f",temp_LRT);
				fprintf(newLOD,"\t%f",temp_LOD);
                        }
                }
		fprintf(newMLE,"\n");
		fprintf(newLRT,"\n");
		fprintf(newLOD,"\n");
		
        }

	fclose(newMLE);
	fclose(newLRT);
	fclose(newLOD);

//重命名
        char oldMLENAME[len + 20];
        char oldLRTNAME[len + 20];
        char oldLODNAME[len + 20];
        strcpy(oldMLENAME,inputPrefix);
        strcpy(oldLRTNAME,inputPrefix);
        strcpy(oldLODNAME,inputPrefix);
        strcat(oldMLENAME,".MLE.old.table.xls");
        strcat(oldLRTNAME,".LRT.old.table.xls");
        strcat(oldLODNAME,".LOD.old.table.xls");

	rename(MLENAME,oldMLENAME);
	rename(LRTNAME,oldLRTNAME);
	rename(LODNAME,oldLODNAME);
	
	rename(newMLENAME,MLENAME);
	rename(newLRTNAME,LRTNAME);
	rename(newLODNAME,LODNAME);

}
	
	
		

