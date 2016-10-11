#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "prepare.h"
#include "order.h"
#include "group.h"

int main(int argc,char *argv[]){
	void mainErrorPrint();
	void allErrorPrint();
	
	int prepare = 0;//等于1时代表需要执行
	int group = 0;
	int order = 0;	
	int all = 0;
	
	if(argc < 2 || argc%2){//没有参数或参数不成对
		mainErrorPrint();
	}	

	if(strcmp(argv[1], "prepare") == 0){
		prepare = 1;
	}else if(strcmp(argv[1], "group") == 0 ){
		group = 1;
		order = 1;
	}else if(strcmp(argv[1], "order") == 0 ){
		order = 1;
	}else if(strcmp(argv[1], "all") == 0){
		prepare = 1;
		group = 1;
		order = 1;
		all = 1;
	}else{
		mainErrorPrint();
	}

	//读取命令行信息
	char *PrepareInputFileName = NULL;
	char GroupInputFileName[1000];
	char *GroupInputFileName_p = NULL;
	char *OrderInputFileName = NULL;
	char *outputFileDir = NULL;
	int grouptype = 1;//群体类型
	char *prefix = NULL;
	int chr = 1;//染色体对数
	int top = 3;//重组率分群输出最优结果的个数
	float accuracy = 0.0001;//重组率计算精度
	int function = 1;//遗传距离函数
	int method = 1;//分群方法参数，1 代表用重组率分群；2 表示用LOD分群；注意CP群体暂时无法计算LOD
	float lod1 = 3;
	float lod2 = 20;//如果采用LOD分群，和joinmap一样
	char *LODInputFileName = NULL;
	char LODInputFileName_char[1000];
	float step = 1;
	

	int x = 0;
	for(x = 2;x < argc;x += 2)
        {
                if(strcmp(argv[x], "-f1") == 0 ){
			PrepareInputFileName = argv[x+1];
		}else if(strcmp(argv[x], "-o") == 0 ){
			outputFileDir = argv[x+1];
		}else if(strcmp(argv[x], "-g") == 0 ){
			grouptype = atoi(argv[x+1]);
		}else if(strcmp(argv[x], "-p") == 0){
			prefix = argv[x+1];
		}else if(strcmp(argv[x], "-f2") == 0){
			strcpy(GroupInputFileName,argv[x+1]);
			GroupInputFileName_p = GroupInputFileName;
		}else if(strcmp(argv[x], "-f") == 0){
			function = atoi(argv[x + 1]);
		}else if(strcmp(argv[x], "-f3") == 0){
			OrderInputFileName = argv[x+1];
		}else if(strcmp(argv[x], "-f4") == 0){
			strcpy(LODInputFileName_char,argv[x+1]);
			LODInputFileName = argv[x+1];
		}else if(strcmp(argv[x], "-c") == 0){
			chr = atoi(argv[x + 1]);
		}else if(strcmp(argv[x], "-t") == 0){
			top = atoi(argv[x + 1]);
		}else if(strcmp(argv[x], "-l1") == 0){
			lod1 = atof(argv[x + 1]);
		}else if(strcmp(argv[x], "-l2") == 0){
			lod2 = atof(argv[x + 1]);
		}else if(strcmp(argv[x], "-s") == 0){
			step = atof(argv[x + 1]);
		}else if(strcmp(argv[x],"-m") == 0){//grouping method
			method = atoi(argv[x + 1]);
		}else if(strcmp(argv[x], "-a") == 0){
			accuracy = atof(argv[x + 1]);
		}else if(strcmp(argv[x], "-h")){
			mainErrorPrint();
		}
        }
	if(outputFileDir == NULL){outputFileDir = "./";}
	if(prefix == NULL){prefix = "prefix";}

	//输出参数信息
	fprintf(stdout,"Weclome to u2map.\n");
	fprintf(stdout,"Version 1.00: released on August 10th, 2015\n");

	//检查参数
	if(all){
		strcpy(LODInputFileName_char,outputFileDir);
		strcat(LODInputFileName_char,prefix);
		char *lod = ".LOD.table.xls";
		strcat(LODInputFileName_char,lod);
		LODInputFileName = LODInputFileName_char;

		if(PrepareInputFileName == NULL){fprintf(stderr,"ERROR : Input file error.Please check \"-f1\" parameter.\n");allErrorPrint();}
		if(grouptype < 1 || grouptype >21){fprintf(stderr,"ERROR : Group type must betweet 1 to 21.\n");allErrorPrint();}
		if(chr < 1){fprintf(stderr,"ERROR : Chromosome's number must >= 1.\n");allErrorPrint();}
		if(top < 1){fprintf(stderr,"ERROR : Group's numbers must >= 1.\n");allErrorPrint();}
		if(method != 1 && method != 2){fprintf(stderr,"ERROR : Grouping method should be 1 or 2;\n");}
		if(function != 1 && function !=2 ){fprintf(stderr,"Maping function must be kosambi or haldane. Please check \"-f\"'s parameter.\n");allErrorPrint();}
	}

	//第一步 prepare
	if(prepare){
		char PrepareOutputFileName[1000];
		strcpy(PrepareOutputFileName,outputFileDir);
		strcat(PrepareOutputFileName,"/");
		strcat(PrepareOutputFileName,prefix);
		fprintf(stdout,"********************\n");
		fprintf(stdout,"Prepare\n");
		fprintf(stdout,"********************\n");
		fprintf(stdout,"Parameters: u2map prepare -f1 %s -o %s -g %d -p %s -a %f\n",PrepareInputFileName,outputFileDir,grouptype,prefix,accuracy);
		step1_prepare(PrepareInputFileName,PrepareOutputFileName,grouptype,accuracy);
		fprintf(stdout,"Prepare finish!\n\n");
		strcpy(GroupInputFileName,PrepareOutputFileName);
		strcat(GroupInputFileName,".MLE.table.xls");
		GroupInputFileName_p = GroupInputFileName;
	}	

	//第二步 分群+排序
	if(group){
		char GroupOutputDir[1000];
		strcpy(GroupOutputDir,outputFileDir);
		strcat(GroupOutputDir,"/group");
		fprintf(stdout,"********************\n");
		fprintf(stdout,"Group\n");
		fprintf(stdout,"********************\n");
		if(method == 1){
			fprintf(stdout,"Parameters: u2map group -f2 %s -m %d -c %d -t %d\n",GroupInputFileName_p,method,chr,top);
			step2_group(GroupInputFileName_p,GroupOutputDir,chr,top,order,function,method,lod1,lod2,LODInputFileName_char,step);
		}else if(method == 2){
			fprintf(stdout,"Parameters: u2map group -f2 %s -f4 %s -m %d -l1 %f -l2 %f -s %f\n",GroupInputFileName_p,LODInputFileName_char,method,lod1,lod2,step);
			step2_group(GroupInputFileName_p,GroupOutputDir,chr,top,order,function,method,lod1,lod2,LODInputFileName_char,step);
		}else{
			
		}
		order = 0;//如果选择的是all ,则group已包含了排序，无需再次排序，
	}

	//第三步 排序+绘图
	if(order){
		step3_order(OrderInputFileName,function);
	}



	return 0;
}

//输入错误提示
void mainErrorPrint(){
	fprintf(stderr,"********************UtopMap********************\n");
	fprintf(stderr,"\nUsage: u2map <command> [option]\n");
	fprintf(stderr,"prepare\t\tPrepare MLE file for grouping or ordering.\n");
	fprintf(stderr,"group\t\tClustering markers into linkage groups.\n");
	fprintf(stderr,"order\t\tOrdering markers in each inkage groups.\n");
	fprintf(stderr,"all\t\tDoing all the above in turn.\n\n");
	exit(EXIT_FAILURE);
}

void allErrorPrint(){
	fprintf(stderr,"Usage:u2map all -f1 mark.input -o . -p prefix -g 1 -c 1 -t 3 -a 0.0001\n");
	fprintf(stderr,"\n\t-f1\t[str] MLE Input file,begin with markNum.\n");
	fprintf(stderr,"\t-o\t[str] Output directory,must exist. [default .]\n");
	fprintf(stderr,"\t-p\t[str] Output file prefix. [default prefix]\n");
	fprintf(stderr,"\t-a\t[float] Accuracy of recombination frequency. [default 0.0001]\n");
	fprintf(stderr,"\t-c\t[int] The number of chromosome. use when method(-m) = 1;[defalut 1]\n");
	fprintf(stderr,"\t-s\t[flost] group step of LOD.[defalut 1]\n");
	fprintf(stderr,"\t-l1\t[float] min LOD for output. use when method(-m) = 2;[defalut 3]\n");
	fprintf(stderr,"\t-l2\t[float] max LOD for output. use when method(-m) = 2;[defalut 20]\n");
	fprintf(stderr,"\t-f\t[int] Mapping function.1 for kosambi, 2 for haldane. [defalut 1]\n");
	fprintf(stderr,"\t-t\t[int] Group numbers to output.[defalut 3]\n");
	fprintf(stderr,"\t-m\t[int] Grouping method.1 for recombine, 2 for LOD.[default 1]\n");
	fprintf(stderr,"\t-h\tPrint this helpful message.\n");
	fprintf(stderr,"\t-g\t[int] Group type, from 1 to 21.\n");
	fprintf(stderr,"\t\t 01:F1DH      02:F1RIL     03:P1BC1RIL \n");
	fprintf(stderr,"\t\t 04:P2BC1RIL  05:P1BC2RIL  06:P2BC2RIL \n");
	fprintf(stderr,"\t\t 07:P1BC1DH   08:P2BC1DH   09:P1BC2DH  \n");
	fprintf(stderr,"\t\t 10:P2BC2DH   11:P1BC1F1   12:P2BC1F1  \n");
	fprintf(stderr,"\t\t 13:F2        14:F3        15:P1BC2F1  \n");
	fprintf(stderr,"\t\t 16:P2BC2F1   17:P1BC1F2   18:P2BC1F2  \n");
	fprintf(stderr,"\t\t 19:P1BC2F2   20:P2BC2F2   21.CP    \n\n");
	exit(EXIT_FAILURE);
}


