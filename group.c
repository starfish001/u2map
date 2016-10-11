/*
20150525 version0.3去掉模拟数据功能,输出到文件
20150611 修正了MinDistList的整合函数，添加新的属性last
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "order.h"
#include "group.h"

//初始化簇，每一个样品为一个簇
void InitCluster(Cluster *clus,int MARK)
{
    int i = 0;
    clus->num = 1;
    Mark *m = (Mark *)malloc(sizeof(Mark));
    m->markno = 0;
    m->nextMark = NULL;
    clus->mark = m;
    clus->lastCluster = NULL;
    clus->similar = NULL;
    Cluster *p1 = clus;    
    for(i=1;i<MARK;i++)
    {
        Cluster *p2 = (Cluster *)malloc(sizeof(Cluster));
	Mark *m = (Mark *)malloc(sizeof(Mark));
	m->markno = i;
	m->nextMark = NULL;
	p2->mark = m;
	p2->num = 1;
	p2->similar = NULL;
	p1->nextCluster = p2;
	p2->lastCluster = p1;
	if(i == MARK - 1){p2->nextCluster = NULL;}else{p1 = p2;}//下一个循环
    }
}

//计算两个簇的距离(最短距离法)
float distMin(Cluster *cluster1,Cluster *cluster2,float *r,int MARK)
{
    int x,y;
    float (*p)[MARK] = (float (*)[MARK])r;
    float min = MAX;
    Mark *m1 = cluster1->mark;
    Mark *m2 = NULL;
    while(m1)
    {
	m2 = cluster2->mark;
	while(m2)
	{
	    x = m1->markno;
	    y = m2->markno;
	    float temp =  *(*(p+x)+y);
	    if(temp < min){min = temp;} 
	    m2 = m2->nextMark;
	}
	m1 = m1->nextMark;
    }
    return min;
}

//计算两个簇的距离(最远距离法)
float distMax(Cluster *cluster1,Cluster *cluster2,float *r,int MARK)
{
    float (*p)[MARK] = (float (*)[MARK])r;
    int x,y;
    float max = -50;
    Mark *m1 = cluster1->mark;
    Mark *m2 = NULL;
    while(m1)
    {
	m2 = cluster2->mark;
        while(m2)
        {
            x = m1->markno;
            y = m2->markno;
	    float temp =  *(*(p+x)+y);
	    if(temp > max){max = temp;} 
            m2 = m2->nextMark;
        }
        m1 = m1->nextMark;
    }
    return max;
}

//计算两个簇的距离(类平均法)
float distAver(Cluster *cluster1,Cluster *cluster2,float *r,int MARK)
{
    float (*p)[MARK] = (float (*)[MARK])r;
    int x,y;
    float N = 0;
    float sum = 0;
    Mark *m1 = cluster1->mark;
    Mark *m2 = NULL;
    while(m1)
    {
	m2 = cluster2->mark;
        while(m2)
        {
            x = m1->markno;
            y = m2->markno;
	    float temp = *(*(p+x)+y);
            sum += pow(temp,2);
	    N++;
            m2 = m2->nextMark;
        }
        m1 = m1->nextMark;
    }
    return sqrt(sum/N);
}

//递归释放Similar 列表内存
void freeSimilar(Similar * si){
    if(si->next){freeSimilar(si->next);}
    si->next =  NULL;
    si->clu = NULL;
    free(si);
    si = NULL;
}

//更新最近簇的距离，并把每一个cluster节点的similar链表补充，（距离最小且相等的cluster节点都在similar上，没有顺序之分），并给簇编号
int updateCluster(Cluster *clus,float *r,int MARK,int method)
{
    Cluster *p1,*p2;
    float mindist,dist;
    Similar *end;
    int i = 0;//簇编号
    p1 = clus;
    int num;

    if(method == 1)
    {
        while(p1)
        {
	    i++;
	    mindist = MAX;
	    p2 = clus;
	    while(p2)
	    {
	        if(p1 != p2)
	        {
	    	    dist = distMin(p1,p2,r,MARK);
	    	    if(dist < mindist)
	    	    {//重置指针数组
	    	        mindist = dist;
	    	    	if(p1->similar){freeSimilar(p1->similar);}
	    	    	Similar *s = (Similar *)malloc(sizeof(Similar));
	    	    	s->clu = p2;
	    	    	s->next = NULL;
	    	    	end = s;//标记末节点
	    	    	p1->similar = s;
	    	    	num = 1;
	    	    }else if(dist == mindist)
	    	    {
	    	    	Similar *s = (Similar *)malloc(sizeof(Similar));
	    	    	s->clu = p2;
	    	    	s->next = NULL;
	    	    	end->next = s;
	    	    	end = s;
	    	    	num++;    
	    	    }
	        }
	        p2 = p2->nextCluster;
	    }
	    p1->similarDistance = mindist;//最近的簇的距离
	    p1->similarNum = num;
	    p1->clusterNo = i;
	    p1 = p1->nextCluster;
        }
    }else
    {
        while(p1)
        {
	    i++;
	    mindist = -50;//LOD分群,记录maxdist 
	    p2 = clus;
	    while(p2)
	    {
	        if(p1 != p2)
	        {
	    	    dist = distMax(p1,p2,r,MARK);
	    	    if(dist > mindist)
	    	    {//重置指针数组
	    	        mindist = dist;
	    	    	if(p1->similar){freeSimilar(p1->similar);}
	    	    	Similar *s = (Similar *)malloc(sizeof(Similar));
	    	    	s->clu = p2;
	    	    	s->next = NULL;
	    	    	end = s;//标记末节点
	    	    	p1->similar = s;
	    	    	num = 1;
	    	    }else if(dist == mindist)
	    	    {
	    	    	Similar *s = (Similar *)malloc(sizeof(Similar));
	    	    	s->clu = p2;
	    	    	s->next = NULL;
	    	    	end->next = s;
	    	    	end = s;
	    	    	num++;    
	    	    }
	        }
	        p2 = p2->nextCluster;
	    }
	    p1->similarDistance = mindist;//最近的簇的距离
	    p1->similarNum = num;
	    p1->clusterNo = i;
	    p1 = p1->nextCluster;
        }

    }
    return i;
 }

//释放MinDistList内存
void freeDistList(MinDistList *minDist){
    minDist->last = NULL;//20150611
    if(minDist->next){freeDistList(minDist->next);}
    if(minDist->sim){freeSimilar(minDist->sim);}//minDist->sim = null
    minDist->next = NULL;
    minDist->sim = NULL;
    free(minDist);
    minDist = NULL;
}

//筛选下一次聚类的簇
float SelectCluster(MinDistList *minDist,Cluster *clus,int method)
{   
    float min;
    Cluster *c = clus;
    MinDistList *minEnd,*node;
    Similar *s1,*s2,*s3;
	    	
    if(method == 1)
    {
	min = MAX;
        while(c)
        {
	    if(c->similarDistance < min)
	    {
	    	freeDistList(minDist);//先清空
	    	minDist = (MinDistList *)malloc(sizeof(MinDistList));
	    	s1 = (Similar *)malloc(sizeof(Similar));
	    	s1->clu = c;
	    	minDist->sim = s1;
	    	minDist->next = NULL;
	    	minDist->last = NULL;//20150611
	    	s2 = s1->clu->similar;
	    	minEnd = minDist;//末节点指针指回头部
	    	while(s2)
	    	{
		    s3 = (Similar *)malloc(sizeof(Similar));//创建新节点
		    s3->clu = s2->clu;
		    s3->next = NULL;
		    s1->next = s3;
		    s1 = s3;
		    s2 = s2->next;
	    	}
	    	s1->next = NULL;//最后一个节点
	    	min = c->similarDistance;//更新最小距离
	    }else if(c->similarDistance == min)
	    {
	    	node = (MinDistList *)malloc(sizeof(MinDistList));
	    	minEnd->next = node;//连接相同距离的链
	    	node->last = minEnd;//20150611
	    	minEnd = node;//指向最后一个链节点
	    	s1 = (Similar *)malloc(sizeof(Similar));
	    	node->sim = s1;
	    	node->next = NULL;
	    	s1->clu = c;
	    	s2 = s1->clu->similar;
            	while(s2)
               {
                    s3 = (Similar *)malloc(sizeof(Similar));//创建新节点
                    s3->clu = s2->clu;
		    s3->next = NULL;
                    s1->next = s3;
                    s1 = s3;
                    s2 = s2->next;
            	}
	    }
	    c = c->nextCluster;
        }
    }
    else
    {
	min = -50;
        while(c)
        {
	    if(c->similarDistance > min)
	    {
	    	freeDistList(minDist);//先清空
	    	minDist = (MinDistList *)malloc(sizeof(MinDistList));
	    	s1 = (Similar *)malloc(sizeof(Similar));
	    	s1->clu = c;
	    	minDist->sim = s1;
	    	minDist->next = NULL;
	    	minDist->last = NULL;//20150611
	    	s2 = s1->clu->similar;
	    	minEnd = minDist;//末节点指针指回头部
	    	while(s2)
	    	{
		    s3 = (Similar *)malloc(sizeof(Similar));//创建新节点
		    s3->clu = s2->clu;
		    s3->next = NULL;
		    s1->next = s3;
		    s1 = s3;
		    s2 = s2->next;
	    	}
	    	s1->next = NULL;//最后一个节点
	    	min = c->similarDistance;//更新最小距离
	    }else if(c->similarDistance == min)
	    {
	    	node = (MinDistList *)malloc(sizeof(MinDistList));
	    	minEnd->next = node;//连接相同距离的链
	    	node->last = minEnd;//20150611
	    	minEnd = node;//指向最后一个链节点
	    	s1 = (Similar *)malloc(sizeof(Similar));
	    	node->sim = s1;
	    	node->next = NULL;
	    	s1->clu = c;
	    	s2 = s1->clu->similar;
            	while(s2)
               {
                    s3 = (Similar *)malloc(sizeof(Similar));//创建新节点
                    s3->clu = s2->clu;
		    s3->next = NULL;
                    s1->next = s3;
                    s1 = s3;
                    s2 = s2->next;
            	}
	    }
	    c = c->nextCluster;
        }

    }
    return min;
}//end

//判断两个Similar链表是否有交集
int isRelated(Similar *p1,Similar *p2)
{
    Similar *t2;
    while(p1)
    {
	t2 = p2;
	while(t2)
	{
	    if(p1->clu->clusterNo == t2->clu->clusterNo){return 1;}
	    t2 = t2->next;
	}
	p1 = p1->next;
    }
    return 0;
}

//合并两个Similar链表，去除重复项；p2的连到p1上
void mergeSimilar(Similar *p1,Similar *p2)
{
    Similar *temp;
    while(p2)
    {
	temp = p1;
	while(temp)
	{
	    if(temp->clu->clusterNo == p2->clu->clusterNo){break;}//重复项忽略
	    if((temp->next) == NULL)//比较完所有p1节点，非重复项；连到p1上
	    {
		Similar *news = (Similar *)malloc(sizeof(Similar));
		news->clu = p2->clu;
		news->next = NULL;
		temp->next = news;
		temp = news;//跳过新增的节点，判断p2的下一个节点
	    }
	    temp = temp->next;
	}
	p2 = p2->next;
    }
}


//整合聚类的簇列表，把有交叉的簇整合成一个簇
void mergeClusterList(MinDistList *minDist)
{
    MinDistList *p1,*p2;
    p1 = minDist;
    while(p1)
    {
        p2 = p1->next;
	while(p2)
	{
	    if(isRelated(p1->sim,p2->sim))//两个similar链表有交集，把P2的sim连到P1的sim，删除P1
	    {
		mergeSimilar(p1->sim,p2->sim);
		/*p1->next = p2->next;//删除p2节点
		p2->next = NULL;//如果不加这条命令，p2->next之后的节点也会被释放*/

		//删除p2节点 修改于20150611
		if(p2->next)//不是最后一个节点
		{	
		    p2->last->next = p2->next;
		    p2->next->last = p2->last;
		    p2->last = NULL;
		    p2->next = NULL;
		}else //最后一个节点
		{
		    p2->last->next = NULL;
		    p2->last = NULL;
		}
		freeDistList(p2);//释放链表内存
		p2 = p1->next;//重新检查
	    }else
	    {
		p2 = p2->next;
	    }
	}
	p1 = p1->next;
    }
}


//根据mergeClusterList得到聚类的簇列表，将Cluster最终整合
//包括MARK LIST整合、节点回收、列表整体链接
Cluster *mergeCluster(MinDistList *minDist,Cluster *clus)
{
    Cluster *clusHead = clus;
    MinDistList *tempMin = minDist;
    while(tempMin)
    {
	Cluster *newCluster = (Cluster *)malloc(sizeof(Cluster));
	//将新节点加到就Cluster的头部
        newCluster->nextCluster = clusHead;
        newCluster->lastCluster =NULL;
	newCluster->mark = NULL;
	newCluster->similar = NULL;
	newCluster->num = 0;
        clusHead->lastCluster = newCluster;
        clusHead = newCluster;
	Similar *tempSimilar = tempMin->sim;
	while(tempSimilar)
	{
	    newCluster->num += tempSimilar->clu->num;
	    if(!newCluster->mark)//第一次复制，直接把指针指向修改
	    {
	        newCluster->mark = tempSimilar->clu->mark;
	    }else//第2次复制及第2次以上复制mark list；将要复制的内容连接到之前的list末尾
	    {
	        Mark *end = newCluster->mark;//标记list末尾
	        while(end->nextMark){end = end->nextMark;}
	        end->nextMark = tempSimilar->clu->mark;//连接到末尾
	    }
	    //释放该节点tempSimilar->clu
	    if(tempSimilar->clu->nextCluster)
	    {//中间节点
	        tempSimilar->clu->lastCluster->nextCluster = tempSimilar->clu->nextCluster;
	        tempSimilar->clu->nextCluster->lastCluster = tempSimilar->clu->lastCluster;
	    }else
	    {//末节点
		tempSimilar->clu->lastCluster->nextCluster = NULL;
	    }
	    tempSimilar->clu->nextCluster = NULL;
	    tempSimilar->clu->lastCluster = NULL;
	    tempSimilar->clu->mark = NULL;
	    tempSimilar->clu->similar = NULL;
	    free(tempSimilar->clu);

	    tempSimilar = tempSimilar->next;
	}
	tempMin = tempMin->next;
    }
    return clusHead;
}

//输出当前的簇
void printCluster(Cluster * clus,char filename1[1000],char filename2[1000],char (*name)[30],float distSubtraction)
{
    FILE *fp1 = fopen(filename1,"w");
    FILE *fp2 = fopen(filename2,"w");
    fprintf(fp1,"distSubtraction/LOD = %.4f\n",distSubtraction);
    fprintf(fp2,"distSubtraction/LOD = %.4f\n",distSubtraction);
    int i = 1;
    while(clus)
    {	
	fprintf(fp1,"cluster%d %d ",i,clus->num);
	fprintf(fp2,"cluster%d %d ",i,clus->num);
	i++;
	Mark *m = clus->mark;
	while(m)
	{
	    fprintf(fp1,"%d,",m->markno);
	    fprintf(fp2,"%s,",name[m->markno]);
	    m = m->nextMark;
	}
        clus = clus->nextCluster;
        fprintf(fp1,"\n");
        fprintf(fp2,"\n");
    }
    fclose(fp1);
    fclose(fp2);
}

/*模拟姐妹染色单体交叉互换
void recombine(float *r,int MARK,int mark[],int num,int start,int end,int INDIV)
{
	float (*p)[MARK] = (float (*)[MARK])r;
	int indiv = 0;
        while(indiv < INDIV)
        {
                int dist = rand() % 100 + 1;
                int i;
                for(i = 0;i < num;i++)
                {
                        if(dist < mark[i])//"dist" was betweent mark[i-1] and mark[i]
                        {
                                if(i==0)break;//"dist" was in front of M0
                                int x,y;
                                for(x = start; x < start + i;x++)
                                {
                                        for(y = start + i;y < end;y++)
                                        {
                                                *(*(p+x)+y) = *(*(p+x)+y) + 1;
                                                *(*(p+y)+x) = *(*(p+y)+x) + 1;
                                        }
                                }
                        break;
                        }
                }
                indiv++;
        }
	
}
*/
/*
//测试模拟生成重组率矩阵
void simulateR(float *r,int MARK,int INDIV,char (*name)[20])
{
    float (*p)[MARK] = (float (*)[MARK])r;
    srand((unsigned)time(NULL));
    int chr1[30] = {2,5,7,10,15,16,21,24,29,30,34,37,42,47,52,54,55,59,62,66,69,70,72,73,83,84,88,90,93,98};
    int chr2[10] = {5,10,22,30,34,55,66,73,84,90};
    int chr3[10] = {5,10,22,30,34,55,66,73,84,90};
    int i = 0;
    int j = 0;
    while(i < MARK)
    {
	j = 0;
	while(j< MARK)
	{
	    *(*(p+i)+j) = 0;
	    j++;
	}
	char temp[9];
	char src[10] = "M";
	sprintf(temp,"%d",i);//将数据格式化的转换成字符串格式
	strcat(src,temp);
	strcpy(*(name+i),src);
	i++;
    }
    recombine(r,MARK,chr1,30,0,30,INDIV);
    recombine(r,MARK,chr2,10,30,40,INDIV);
    recombine(r,MARK,chr3,10,40,50,INDIV);
    int x,y;
    printf("markNum=%d\nMARK\t",MARK);
    for(x=0;x<MARK;x++){printf("M%d\t",x);}
    printf("\n");
    for(x = 0;x<MARK;x++)
    {
	printf("M%d\t",x);
	for(y = 0; y<MARK;y++)
	{
	    float temp = *(*(p+x)+y);
	    if(temp > 0)
	    {
	        float temp = *(*(p+x)+y);
	        temp = temp / INDIV / 2;
	        //temp = 25 * log((1 + 2 * temp) / (1 - 2 * temp));
	        *(*(p+x)+y) = temp;
	    }else
	    {
	        *(*(p+x)+y) = 0.75;
	    }
	    printf("%.2f\t",*(*(p+x)+y));
	}
	printf("\n");
    }
}
*/
//通过读文件获取重组率的矩阵表
void readTable(char file[],float *r,char (*name)[30])
{
    FILE *fp;
    fp = fopen(file,"r");
    int markNum = 0;//MARK的数量

    fscanf(fp,"markNum=%d",&markNum);

    float (*p)[markNum] = (float (*)[markNum])r;

    int i,j;
    for(i = 0;i < markNum;i++)
    {
        fscanf(fp,"%s",*(name+i));//读取MARK name
        for(j = i;j < markNum;j++)
        {
	    if(i == j)
	    {
		*(*(p+i)+j) = 0.75;
	    }else
	    {
	        float temp;
                fscanf(fp,"%f",&temp);
	        *(*(p+i)+j) = temp;
	        *(*(p+j)+i) = temp;
	    }
        }
    }
    fclose(fp);
}


//读取LOD文件
void readLODTable(char file[],float *r)
{
    FILE *fp;
    fp = fopen(file,"r");
    int markNum = 0;//MARK的数量

    fscanf(fp,"markNum=%d",&markNum);

    float (*p)[markNum] = (float (*)[markNum])r;

    int i,j;
    char name[30];
    for(i = 0;i < markNum;i++)
    {
        fscanf(fp,"%s",name);//读取MARK name
        for(j = i;j < markNum;j++)
        {
	    if(i == j)
	    {
		*(*(p+i)+j) = 0;
	    }else
	    {
	        float temp;
                fscanf(fp,"%f",&temp);
		//printf("%f\n",temp);
		if(isnan(temp))
		{
			//printf("not a num.\n");
			temp=0;
		}
	        *(*(p+i)+j) = temp;
	        *(*(p+j)+i) = temp;
	    }
        }
    }
}

//输入当前聚类分群结果到文件中
void printToFile(char *inputfile,char *outputDir,float *dist,char (*name)[30], int MARK,float index,int order,int function,int method)
{
	float (*r)[MARK] = (float (*)[MARK])dist;	

	//新建文件夹
	if(strlen(outputDir) > 1000){fprintf(stderr,"ERROR: the length of file name is longer than 1000.\n");exit(EXIT_FAILURE);}
	char command[1036] = "mkdir ";
	char dir[1000];
	char n[4];
	sprintf(n,"%.2f",index);//将数据格式化的转换成字符串格式
	strcpy(dir,outputDir);
	if(method == 1)
	{
		strcat(dir,"/TOP");
	}else
	{
		strcat(dir,"/LOD_");
	}
	strcat(dir,n);
	strcat(dir,".Groups");
	strcat(command,dir);
	system(command);

	FILE *in = NULL;
	in =fopen(inputfile,"r");
	if(in == NULL){fprintf(stderr,"%s input error.\n",inputfile);exit(EXIT_FAILURE);}
	float distance;
	fscanf(in,"distSubtraction/LOD = %f",&distance);
	while(!feof(in))
	{
		char cname[13];
		int num;
		//读取文件
		fscanf(in,"%s ",cname);
		fscanf(in,"%d ",&num);
		int MarkIndex[num];
		int i = 0;
		while(i < num){
			fscanf(in,"%d,\n",&MarkIndex[i]);
			i++;
		}
		//输出文件
		char file[1120];
		strcpy(file,dir);
		strcat(file,"/");
		strcat(file,cname);
		strcat(file,".MLE.xls");

		FILE *fp;
		fp = fopen(file,"w");
		i = 0;
		fprintf(fp,"markNum=%d",num);//该簇包含Mark数
		while(i < num)
		{
			int x = MarkIndex[i];
			fprintf(fp,"\n%s",name[x]);
			int j = i + 1;
			while(j < num)
			{
				int y =  MarkIndex[j];
				fprintf(fp,"\t%f",*(*(r+x)+y));
				j++;
			}
			i++;
		}
		fclose(fp);
		if(order){step3_order(file,function);}//排序
	}
	fclose(in);
	
}

//错误输出
void GroupErrorPrint()
{
	fprintf(stderr,"Usage1:u2map group -f2 MLE.xls -c 10 -t 3 -m 1\n");
	fprintf(stderr,"\n\t-f2\t[str] MLE file.\n");
	fprintf(stderr,"\t-m\t[int] Group method. 1 for MLE, 2 for LOD.[defalut 1]\n\n");
	fprintf(stderr,"\t-c\t[int] The number of chromosome. set when method(-m) = 1;[defalut 1]\n");
	fprintf(stderr,"\t-t\t[int] Group numbers to output. set when method(-m) = 1;[defalut 3]\n\n");
	fprintf(stderr,"\n\t-f4\t[str] LOD Input file,begin with markNum. set when method(-m) = 2;\n");
	fprintf(stderr,"\t-l1\t[float] min LOD for output. set when method(-m) = 2;[defalut 3]\n");
	fprintf(stderr,"\t-l2\t[float] max LOD for output. set when method(-m) = 2;[defalut 20]\n");
	fprintf(stderr,"\t-s\t[flost] group step of LOD. set when method(-m) = 2;[defalut 1]\n\n");

	exit(EXIT_FAILURE);
}

int step2_group(char *file,char *output,int CHR,int TOP_Num,int order,int function,int method,float lod1,float lod2,char *lodfile,float step)
{
	int MARK = 0;//mark 数量
	//int CHR = 0;//染色体数量，可有参数设定；默认为1
	//int TOP_Num = 0;//分群过程中，在分群数不小于染色体数量的前提下，输出跨度最大的3次聚类，默认是3
	//char *output = NULL;//输出文件的前缀

	if(file == NULL){GroupErrorPrint();}
	FILE *fp;
	if((fp = fopen(file,"r")) != NULL){
		fscanf(fp,"markNum=%d",&MARK);
		fclose(fp);
	}else{
		fprintf(stderr,"ERROR: can not open %sfile.",file);
		exit(EXIT_FAILURE);	
	}
	//验证、设定默认参数
	if(CHR == 0){CHR = 1;}
	if(MARK == 0){fprintf(stderr,"ERROR: MARK NUM = 0;pleasa check input file.\n");GroupErrorPrint();}
	if(TOP_Num == 0){TOP_Num = 3;}
	if(!output){fprintf(stderr,"ERROR: output file prefix wrong.\n");GroupErrorPrint();}
	if(method != 1 && method != 2){fprintf(stderr,"ERROR: method must be 1 or 2.");GroupErrorPrint();}
	if(lodfile == NULL){fprintf(stderr,"ERROR,f4 can not be null.\n");GroupErrorPrint();}

	//开始读MLE文件
	float dist[MARK][MARK];
	char markName[MARK][30];//存储标记的名称
	float *r = (float *)dist;
	readTable(file,r,markName);

	//创建目录
	char common[1000] = "mkdir ";
	strcat(common,output);
	system(common);
	
	//MLE分群
	if(method == 1)
	{
		//聚类过程中，先输出相应Top3的聚类表，边聚类，边输出。
		char clusterFile[TOP_Num][1000];//每个聚类结果的文件名，一共TOP_Num个
		char clusterNameFile[TOP_Num][1000];//每个聚类结果的文件名，一共TOP_Num个
		float Topdist[TOP_Num];//保存最远距离的几个距离变量
		float lastDist;//上一次的聚类距离

		int fi = 0;
		while(fi < TOP_Num){
			Topdist[fi] = 0.0;
			char name1[1000];
			char name2[1000];
			strcpy(name1,output);
			strcat(name1,"/TOP");
			char n[4];
			sprintf(n,"%d",fi+1);//将数据格式化的转换成字符串格式
			strcat(name1,n);
			strcpy(name2,name1);
			strcat(name1,".cluster.index.txt");
			strcat(name2,".cluster.name.txt");
			strcpy(clusterFile[fi],name1);
			strcpy(clusterNameFile[fi],name2);
			FILE *fp1;
			FILE *fp2;
			fp1 = fopen(name1,"w");
			fp2 = fopen(name2,"w");
			fclose(fp1);
			fclose(fp2);
			fi++;
		}//文件生成
		
		//过程输出信息
		FILE *infor;
		char logName[1000];
		strcpy(logName,output);
		strcat(logName,"/Group.log");
		infor = fopen(logName,"w");
	
		Cluster *ClusHead = (Cluster *)malloc(sizeof(Cluster));
		MinDistList *minDist = (MinDistList *)malloc(sizeof(MinDistList));
		ClusterDist *distListHead = (ClusterDist *)malloc(sizeof(ClusterDist));
		ClusterDist *distListEnd = NULL;
		distListHead->next = NULL;
		distListHead->cluNum = MARK;
		minDist->sim = NULL;
		minDist->next = NULL;
		minDist->last = NULL;//20150611
		
		//只需一次聚类，但是分群过程不断刷新输出;
		InitCluster(ClusHead,MARK);
		updateCluster(ClusHead,r,MARK,method);
		//开始聚类，输出聚类的结果范围是最终聚类的簇数在染色体条数和2倍染色体条数之间（CHR ~ 2*CHR）
		while(ClusHead->nextCluster){//最后成为一簇,next为空
			int clusterNUM;//当前聚类的簇数
			if(!distListEnd)//第一次聚类
			{
				distListHead->dist = SelectCluster(minDist,ClusHead,method);
				lastDist = distListHead->dist;
				mergeClusterList(minDist);
				ClusHead = mergeCluster(minDist,ClusHead);
				clusterNUM = updateCluster(ClusHead,r,MARK,method);
				distListHead->cluNum = clusterNUM;
				distListEnd = distListHead;
			}else//第n次聚类
			{
				float thisDist;
				ClusterDist *temp = (ClusterDist *)malloc(sizeof(ClusterDist));
				temp->dist = SelectCluster(minDist,ClusHead,method);
				thisDist = temp->dist;
				mergeClusterList(minDist);
				ClusHead = mergeCluster(minDist,ClusHead);
				clusterNUM = updateCluster(ClusHead,r,MARK,method);
				temp->cluNum = clusterNUM;
				temp->next =NULL;
				distListEnd->next = temp;
				distListEnd = temp;
	
				int i = TOP_Num - 1;
				float distSubtraction = thisDist - lastDist;
				int flag;
				if(clusterNUM >=  CHR && distSubtraction > Topdist[i]){
					i--;
					while(i >= 0){
						if(distSubtraction > Topdist[i]){
							rename(clusterFile[i],clusterFile[i+1]);
							rename(clusterNameFile[i],clusterNameFile[i+1]);
							Topdist[i+1] = Topdist[i];
							if(i == 0){flag = i;Topdist[i] = distSubtraction;}
						}else{
							flag = i + 1;
							Topdist[i + 1] = distSubtraction;
							break;
						}
						i--;
					}
					printCluster(ClusHead,clusterFile[flag],clusterNameFile[flag],markName,Topdist[flag]);
				}
				lastDist = thisDist;
			}
		}
	
		//输出详细聚类分群的结果
		fprintf(stdout,"Group finish!\n\n");
		int i;
		for(i = 0;i < TOP_Num;i++){
			printToFile(clusterFile[i],output,r,markName,MARK,i+1,order,function,method);
			remove(clusterFile[i]);
		}
	
		//输出聚类过程信息
		ClusterDist *distp = distListHead;
		int times = 1;
		while(distp)
		{
			fprintf(infor,"cluster_times: %d\tDist: %.4f\tPercentage: %.4f%%\tcluster_num: %d\n",times,distp->dist,(distp->dist) / (distListEnd->dist) * 100,distp->cluNum);
			distp = distp->next;
			times++;
		}
		
		//printCluster(ClusHead,markName);
		
		fclose(infor);
	}else //LOD分群
	{
		//LOD 最多输出30个文件
		char clusterFile[100][1000];//每个聚类结果的文件名，一共TOP_Num个
		char clusterNameFile[100][1000];//每个聚类结果的文件名，一共TOP_Num个
		float lod_no[100];//每个群的lod值
		/*
		int fi = lod2;
		int filenum = 0;
		while(fi >=  lod1){
			char name1[1000];
			char name2[1000];
			strcpy(name1,output);
			strcat(name1,"/LOD");
			char n[4];
			sprintf(n,"%d",fi+1);//将数据格式化的转换成字符串格式
			strcat(name1,n);
			strcpy(name2,name1);
			strcat(name1,".cluster.index.txt");
			strcat(name2,".cluster.name.txt");
			strcpy(clusterFile[filenum],name1);
			strcpy(clusterNameFile[filenum],name2);
			//FILE *fp1;
			//FILE *fp2;
			//fp1 = fopen(name1,"w");
			//fp2 = fopen(name2,"w");
			//fclose(fp1);
			//fclose(fp2);
			fi--;
			filenum++;
		}//文件生成
		*/
		//开始读MLE文件
		float lod_table[MARK][MARK];
		float *lod_p = (float *)lod_table;
		readLODTable(lodfile,lod_p);

		Cluster *ClusHead = (Cluster *)malloc(sizeof(Cluster));
		MinDistList *maxDist = (MinDistList *)malloc(sizeof(MinDistList));
		maxDist->sim = NULL;
		maxDist->next = NULL;
		maxDist->last = NULL;//20150611
		InitCluster(ClusHead,MARK);
		int clusterNum = updateCluster(ClusHead,lod_p,MARK,method);
		float tempLOD = SelectCluster(maxDist,ClusHead,method);

		int flag = 0;
		while(clusterNum > 1 && ((tempLOD=SelectCluster(maxDist,ClusHead,method)) >= lod1))
		{
			mergeClusterList(maxDist);
			ClusHead = mergeCluster(maxDist,ClusHead);
			if(tempLOD <= lod2 )
			{
				char name1[1000];
				char name2[1000];
				strcpy(name1,output);
				strcat(name1,"/LOD_");
				char n[4];
				
				sprintf(n,"%.2f",tempLOD);//将数据格式化的转换成字符串格式
				strcat(name1,n);
				strcpy(name2,name1);
				strcat(name1,".cluster.index.txt");
				strcat(name2,".cluster.name.txt");
				strcpy(clusterFile[flag],name1);
				strcpy(clusterNameFile[flag],name2);
				printCluster(ClusHead,clusterFile[flag],clusterNameFile[flag],markName,tempLOD);
				//lod2 = (int)tempLOD - step;
				lod2 = tempLOD - step;
				lod_no[flag] = tempLOD;
				flag++;
			}
			clusterNum = updateCluster(ClusHead,lod_p,MARK,method);
		}

		//输出详细聚类分群的结果
		fprintf(stdout,"Group finish!\n\n");
		int i;
		for(i = 0;i < flag;i++){
			printToFile(clusterFile[i],output,r,markName,MARK,lod_no[i],order,function,method);
			remove(clusterFile[i]);
		}
	}
	return 0;
}
