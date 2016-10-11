/*
20150701 排序时，先将重组率为0的marks整合成一个bin，整合之后的重组率按最小的替换

20150609 修改了2-opt和重组算法，当新重组的遗传距离为0时，则放弃该重组。因为当实际数据中很多重组率为0时，有些mark对会多次重组，进入死循环
	 迭代2-opt和重组算法，知道都不能缩小权重
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "order.h"

struct markNode *pathRoot = NULL;
int pathLength = 0;
int UselessMark = 0;//distance = 0的marks数

void MST(float *dist,int markNUM,MARKNode *root)
{
    float (*r)[markNUM] = (float (*) [markNUM])dist;
    int count = 1;//用于记录已建树的节点数
    //r[][]的第一行的MARK作为最小生成树的树根
    root->markIndex = 0;
    root->Level = 0;
    root->parent = NULL;
    root->leftChild = NULL;
    root->rightSibling = NULL;
    
    //listHead用于保存建树过程中，已建树的节点
    TreeList *listHead = (TreeList *)malloc(sizeof(TreeList));
    listHead->markIndex = 0;
    listHead->mstNode = root;
    listHead->next = NULL;
    TreeList *listEnd = listHead;

    //开始构建最小生成树，用左孩子右兄弟方法表示树
    extern int UselessMark;
    int NUM = markNUM - UselessMark;
    //printf("UselessMark=%d,NUM=%d\n",UselessMark,NUM);
    while(count < NUM)
    {
	TreeList *list = listHead;
	TreeList *minNode = listHead;//用于标记listHead选中的节点，从而根据minNode->mstNode定位到最小生成树
	int index1,index2;
	float minDist = MAXRLE;

	//找到下一次建树的节点
	while(list)
	{
	    int x,y;
	    x = list->markIndex;
	    for(y = 0;y < markNUM; y++)
	    {
		if( *(*(r+x)+y) < minDist)
		{
		    minDist = *(*(r+x)+y);
		    minNode = list;
		    index1 = x;
		    index2 = y;
		}
	    }
	    list = list->next;
	}

	//为了确保选取的index2还没建树，当选取了下一个index2值时，需将index2和TreeList的所有节点的R值变成1000，确保下次不会被选取到
	list = listHead;
	while(list)
	{
	    *(*(r + list->markIndex) + index2) = MAXRLE;
	    *(*(r + index2) + list->markIndex) = MAXRLE;
	    list = list->next;
	}

	//在MST上插入节点
	MARKNode *parent = minNode->mstNode;
	MARKNode *newMST = (MARKNode *)malloc(sizeof(MARKNode));
	newMST->markIndex = index2;
	newMST->Level = 0;
	newMST->parent = parent;
	newMST->leftChild = NULL;
	newMST->rightSibling =NULL;
	if(parent->leftChild)
	{//当已将有左孩子是，新增节点作为左孩子的右兄弟
		MARKNode *right = parent->leftChild;
		while(right->rightSibling){right = right->rightSibling;}
		right->rightSibling = newMST;
	}else
	{//否则直接成为左孩子
		parent->leftChild = newMST;
	}

	//将节点加入链表
	listEnd->next = (TreeList *)malloc(sizeof(TreeList));
	listEnd->next->markIndex = index2;
	listEnd->next->mstNode = newMST;
	listEnd->next->next = NULL;
	listEnd = listEnd->next;

	count++;
    }
}

float MLE2MapingDistance(float r,int function){
	if(function == 1){//Kosambi作图函数
		return 25 * log((1 + 2 * r) / (1 - 2 * r));
	}else if(function == 2){//Haldane作图函数
		return log(1 - 2 * r) * (-50);
	}else{
		fprintf(stderr,"Maping function must be kosambi or haldane. Please check \"-f\"'s parameter.\n");exit(EXIT_FAILURE);
	}

}

void readGroupTable(char file[],float *r,char (*name)[30],int function)
{
    FILE *fp;
    int markNum = 0;//MARK的数量
    if((fp = fopen(file,"r"))!= NULL){
        fscanf(fp,"markNum=%d",&markNum);
    }else{
	fprintf(stderr,"can not open %s file.\n",file);
	exit(EXIT_FAILURE);
    }

    float (*p)[markNum] = (float (*)[markNum])r;

    int i,j;
    for(i = 0;i < markNum;i++)
    {
	fscanf(fp,"%s",*(name+i));//读取MARK name
        for(j = i;j < markNum;j++)
        {
            if(i == j)
            {
                *(*(p+i)+j) = MAXRLE;
            }else
            {
                float temp;
                fscanf(fp,"%f",&temp);
		if(temp >= 0.5){temp = 1000;}else{temp = MLE2MapingDistance(temp,function);} //Kosambi作图函数更好地表示重组率
                *(*(p+i)+j) = temp;
                *(*(p+j)+i) = temp;
            }

        }
    }
    fclose(fp);
}


//整和重组率矩阵，将重组率为0的marks整合在一起。保留其中一个，将其他的mark重组率替换成2000,保留的是index小的marks
BinNode * mergeTable(float *dist,int markNum){
	float (*r)[markNum] = (float (*) [markNum])dist;
	extern int UselessMark;
	UselessMark = 0;
	int ii,jj;
	BinNode * Bhead = NULL;
	BinNode * Bend = NULL;
	for(ii=0;ii<markNum;ii++){
		for(jj=ii+1;jj<markNum;jj++){
			if(*(*(r+ii)+jj) == 0){//第ii行和第jj行的mark整合，以及第ii列和第jj列的marks整合，ii行、列用小的替换，jj行列替换成UNLESS
				//先判断ii是否已经包含在链表中，如果是，则将yy加入otherMarks，否则新建一个nextBin
				UselessMark++;
				BinNode * temp = Bhead;
				while(temp){
					if(temp->keepMark == ii){break;}
					Bend = temp;
					temp = temp->nextBin;
				}

				if(temp){//temp不为空，说明ii已经包含在链表中
					Path *path = (Path *)malloc(sizeof(Path));
					path->markIndex = jj;
					path->next = NULL;
					path->last = NULL;

					Path *end = temp->otherMarks;
					while(end->next){end = end->next;}
					end->next = path;
					path->last=end;
				}else{//temp为空，说明未找到ii，新建一个BinNode节点
					Path *path = (Path *)malloc(sizeof(Path));
					path->markIndex = jj;
					path->next = NULL;
					path->last = NULL;
			
					BinNode *newbin = (BinNode *)malloc(sizeof(BinNode));
					newbin->keepMark = ii;
					newbin->otherMarks = path;
					newbin->nextBin=NULL;
					if(!Bhead){//第一个节点
						Bhead = newbin;
					}else{
						Bend->nextBin = newbin;
					}
				}

				int xx = 0;
				while(xx < markNum){
					if( ii != xx && *(*(r+ii)+xx) > *(*(r+jj)+xx)){ *(*(r+ii)+xx) = *(*(r+jj)+xx);}//ii行替换成小的
					if( ii != xx && *(*(r+xx)+ii) > *(*(r+xx)+jj)){ *(*(r+xx)+ii) = *(*(r+xx)+jj);}//ii列替换成小的
					*(*(r+jj)+xx) = UNLESS;//jj行替换成UNLESS
					*(*(r+xx)+jj) = UNLESS;//jj列替换成UNLESS
					xx++;
				}
			}
		}
	}
	return Bhead;
	
}

//更新最小生成树的LEVEL
int modifyLevel(MARKNode *node)
{
    if(node->leftChild)
    {
	int max = 0;
	MARKNode *temp = node->leftChild;
	int level = modifyLevel(temp);
	if(level > max) max = level;
        while(temp->rightSibling)
	{
	    level = modifyLevel(temp->rightSibling);
	    if(level > max) max = level;
	    temp = temp->rightSibling;
	}
	node->Level = max + 1;
	return max+1;
    }else
    {
	node->Level = 1;
	return 1;
    }
}

//输出最小生成树,测试用
int printMST(MARKNode *node,char (*name)[30])
{
    printf("%s:",*(name + node->markIndex));
    if(node->parent)printf("p %s\t",*(name + node->parent->markIndex));
    if(node->leftChild)printf("l %s\t",*(name + node->leftChild->markIndex));
    if(node->rightSibling)printf("r %s\t",*(name + node->rightSibling->markIndex));
    printf("level:%d",node->Level);
    printf("\n");
    MARKNode *temp = node->leftChild;
    while(temp)
    {
	printMST(temp,name);
        temp = temp->rightSibling;
    }
    return 1;
}


//返回最大的孩子指针
MARKNode * maxChild(MARKNode *node)
{
    MARKNode *child = node->leftChild;
    if(child)
    {
	int max = child->Level;
	MARKNode *maxChild = child;
	while(child->rightSibling)
	{
	    child = child->rightSibling;
	    if(child->Level > max)
	    {
		max = child->Level;
		maxChild = child;
	    }
	}
	return maxChild;
    }else
    {//没有孩子，返回空
	return NULL;
    }
}

//从MST上选出最长路径的Middle节点
void getLongestPathRoot(MARKNode *root)
{
    extern struct markNode *pathRoot;
    pathRoot = NULL;
    extern int pathLength;
    pathLength = 0;
    int max1 = 0;//最大子树层数
    int max2 = 0;//次大子树层数
    int max;//最大路径，等于max1+max2+1

    MARKNode *temp = root->leftChild;
    while(temp)
    {
	if(max1 == 0)
	{
	    max1 = temp->Level;//第一个孩子
	}
	else
	{
	    int level = temp->Level;
	    if(level > max2 && level < max1){max2 = level;}
            else if(level > max1){max2 = max1;max1 = level;}
	}
	getLongestPathRoot(temp);//递归遍历
        temp = temp->rightSibling;
		
    }
    max = max1 + max2 + 1;
    if(max > pathLength)
    {
        pathRoot = root;
        pathLength = max;
    }
}

//获取未排序的mark
Path * getNoOrder(Path *orderPath,int markNum,float *dist)
{
    float (*r)[markNum] = (float (*) [markNum])dist;
    int i,flag;
    Path *head =NULL;
    Path *end =NULL;
    for(i = 0;i<markNum;i++)
    {
	if( *(*(r+i)+0) == UNLESS ){continue;}//该mark已被整合
        flag = 0;//如果flag为1，表示编号为i的标记已经在orderPath中
        Path *temp = orderPath;
        while(temp)
        {
            if(temp->markIndex == i){flag = 1;break;}
	    temp = temp->next;
        }
        if(flag){continue;}//下一个i
        else
        {
            Path *temp = (Path *)malloc(sizeof(Path));
            temp->markIndex = i;
            temp->next = NULL;
	    if(!head)//第一个标记作为头结点
	    {
	        temp->last = NULL;
		head = temp;
		end = temp;
	    }
            else
	    {
		temp->last = end;
		end->next = temp;
        	end = temp;
	    }	
        }
    }
    return head;
}

/*获取MST的最长路径
root位最长路径的根节点,从root找到最长和次长的子节点串，root在中间，将两个字串连成Path
*/
Path * selectLongestPath(MARKNode *root)
{
    Path *orderPath = (Path *)malloc(sizeof(Path));
    Path *end = (Path *)malloc(sizeof(Path));
    end->markIndex = root->markIndex;
    end->next = NULL;
    end->last = NULL;
    MARKNode *max = maxChild(root);
    Path *head = end;
    max->Level = 0;//此举为了下次选root最大child时不再选到这个节点
    while(max)
    {
	if(max->leftChild)
	{
    	    Path *temp = (Path *)malloc(sizeof(Path));
	    temp->markIndex = max->markIndex;
	    temp->next = head;
	    head->last = temp;
	    temp->last = NULL;
	    head = temp;
	    max = maxChild(max);
	}else
	{//最后一个节点，此节点作为path的头结点,无需再开辟内存
	    orderPath->markIndex = max->markIndex;
	    orderPath->next = head;
	    head->last = orderPath;
	    orderPath->last = NULL;
	    break;
	}
    }
    max = maxChild(root);//选取root的次大链表
    if(max->Level == 0)return orderPath;//finish! root节点只有一个孩子,该孩子的level在上面代码已将其置0
    while(max)
    {
        Path *temp = (Path *)malloc(sizeof(Path));
        temp->markIndex = max->markIndex;
        temp->next = NULL;
	temp->last = end;
        end->next = temp;
	end = temp;
        max = maxChild(max);
    }
    return orderPath;
}

//判断是否为整合的节点,如何是则输出markname和删除该节点
BinNode * isBin(int index,BinNode *head,float distance,char (*name)[30],char *OrderOutputXLS){
	FILE *xls;
	xls = fopen(OrderOutputXLS,"a");
	BinNode * p1 = head;
	if(p1->keepMark == index){
		head = p1->nextBin;
		p1->nextBin = NULL;
		Path * n1 = p1->otherMarks;
		while(n1){
			fprintf(xls,"%s\t%f\t%f\n",*(name + n1->markIndex),0.0,distance);
			n1 = n1->next;
		}
		free(p1);
		p1 = NULL;
		fclose(xls);
		return head;
	}else{
		BinNode *p2 = p1->nextBin;
		while(p2){
			if(p2->keepMark == index){
				Path * n1 = p2->otherMarks;
				while(n1){
					fprintf(xls,"%s\t%f\t%f\n",*(name + n1->markIndex),0.0,distance);
					n1 = n1->next;
				}
				p1->nextBin = p2->nextBin;
				p2->nextBin = NULL;
				free(p2);
				p2 = NULL;
				break;
			}else{
				p1 = p2;
				p2 = p2->nextBin;
			}
		}
		fclose(xls);
		return head;
	}
}

//输出path的内容
void printPath(Path *path,char (*name)[30],float *dist,int markNum,BinNode *useless,char *OrderOutputXLS)
{
    FILE *xls;
    xls = fopen(OrderOutputXLS,"w");
    float (*r)[markNum] = (float (*) [markNum])dist;
    float total = 0;
    while(path)
    {
	int x = path->markIndex;
	if(path->last)//非第一个节点
	{
	    int y = path->last->markIndex;
	    float distance = *(*(r+x)+y);
	    total += distance;
	    fprintf(xls,"%s\t%f\t%f\n",*(name + x),distance,total);
	    if(useless){
		    fclose(xls);
		    useless =  isBin(x,useless,total,name,OrderOutputXLS);
		    xls = fopen(OrderOutputXLS,"a");
	    }
	}else//第一个节点
	{
	    fprintf(xls,"%s\t%f\t%f\n",*(name + x),total,total);
	    if(useless){
		    fclose(xls);
		    useless =  isBin(x,useless,total,name,OrderOutputXLS);
		    xls = fopen(OrderOutputXLS,"a");
		}
	}
	path =  path->next;
    }
    fclose(xls);

}

/*
将noOrderPath的mark根据最小权重和整合到orderPath
如果noOrderPath为空则无需执行此步骤
先找到中间插入位点，然后比较链头插入和链尾插入，确定最终的插入位点
*/
Path * mergePath(Path *orderHead,Path *noOrderHead,float *dist,int markNum)
{
    float (*r)[markNum] = (float (*) [markNum])dist;
    Path *noOrderTemp = noOrderHead;
    Path *Head = orderHead;
    Path *End = NULL;
    while(noOrderTemp)
    {
	Path *p1 = Head;
	Path *p2 = Head->next;
	Path *minNode = NULL;
	
	float min = 2000;//min = r1 + r2 - r3;r1和r2去最大值1000，r3取最小值0
	while(p2)
	{
	    float r1,r2,r3;
	    r1 = *(*(r + p1->markIndex) + noOrderTemp->markIndex);
	    r2 = *(*(r + p2->markIndex) + noOrderTemp->markIndex);
	    r3 = *(*(r + p1->markIndex) + p2->markIndex);
	    float temp = r1 + r2 - r3;
	    if(temp < min)
	    {
		min = temp;
		minNode = p1;//minNode的后面插入新节点
	    }
	    p1 = p1->next;
	    p2 = p2->next;
	}
	End = p1;//链表末尾
	
	//插入新节点
	float headValue = *(*(r + Head->markIndex) + noOrderTemp->markIndex);//新节点与链头的权重
	float endValue = *(*(r + End->markIndex) + noOrderTemp->markIndex);//新节点与链尾的权重
	Path *new = (Path *)malloc(sizeof(Path));
	new->markIndex = noOrderTemp->markIndex;
	if(min < headValue && min < endValue)//从中间插入
	{
	    new->last = minNode;
	    new->next = minNode->next;
	    new->next->last = new;
	    new->last->next = new;
	}else if(headValue <= min && headValue <= endValue)//从头部插入
	{
	    Head->last = new;
	    new->next = Head;
	    new->last = NULL;
	    Head = new;
	}else if(endValue <= min && endValue <= headValue)//从尾部插入
	{
	    new->next = NULL;
	    new->last = End;
	    End->next = new;
	    End = new;
	}
    noOrderTemp = noOrderTemp->next;
    }
    return Head;
}

//倒序p1到p2之间的节点,p1是原表头，p2是原表尾
void invertedOrder(Path *p1,Path *p2)
{
    Path *p = p1->next;
    while(p != p2)
    {
	Path *tmp = p->last;
        p->last = p->next;
        p->next = tmp;
        p = p->last;
    }
    //改变表头和表尾
    p2->next = p2->last;
    p2->last = NULL;
    p1->last = p1->next;
    p1->next = NULL;
}

/*
优化路径
选取i和j两个节点时，将Pi连到Pj；将Pi+1连到Pj+1；从Pi+1到Pj的节点顺序变成倒序
*/
Path * changePath(Path *p1,Path *p2,Path *p3,Path *p4,Path *p5,Path *p6,int type)
{
    Path *head;

    //切成三部分
    p2->next = NULL;
    p3->last = NULL;
    p4->next = NULL;
    p5->last = NULL;
    
    //根据type倒序后连接三部分
    switch(type)
    {
	case  2 : invertedOrder(p3,p4);							head = p1;p2->next = p4;p4->last = p2;p3->next = p5;p5->last = p3;break;
	case  3 : invertedOrder(p5,p6);							head = p1;p2->next = p3;p3->last = p2;p4->next = p6;p6->last = p4;break;
	case  4 : invertedOrder(p3,p4);invertedOrder(p5,p6);				head = p1;p2->next = p4;p4->last = p2;p3->next = p6;p6->last = p3;break;

	case  5 : invertedOrder(p1,p2);							head = p2;p1->next = p3;p3->last = p1;p4->next = p5;p5->last = p4;break;
	case  6 : invertedOrder(p1,p2);invertedOrder(p3,p4);				head = p2;p1->next = p4;p4->last = p1;p3->next = p5;p5->last = p3;break;
	case  7 : invertedOrder(p1,p2);invertedOrder(p5,p6);				head = p2;p1->next = p3;p3->last = p1;p4->next = p6;p6->last = p4;break;
	case  8 : invertedOrder(p1,p2);invertedOrder(p3,p4);invertedOrder(p5,p6);	head = p2;p1->next = p4;p4->last = p1;p3->next = p6;p6->last = p3;break;

	case  9 : 									head = p3;p4->next = p1;p1->last = p4;p2->next = p5;p5->last = p2;break;
	case 10 : invertedOrder(p1,p2);							head = p3;p4->next = p2;p2->last = p4;p1->next = p5;p5->last = p1;break;
	case 11 : invertedOrder(p5,p6);							head = p3;p4->next = p1;p1->last = p4;p2->next = p6;p6->last = p2;break;
	case 12 : invertedOrder(p1,p2);invertedOrder(p5,p6);				head = p3;p4->next = p2;p2->last = p4;p1->next = p6;p6->last = p1;break;

	case 13 : invertedOrder(p3,p4);							head = p4;p3->next = p1;p1->last = p3;p2->next = p5;p5->last = p2;break;
	case 14 : invertedOrder(p1,p2);invertedOrder(p3,p4);				head = p4;p3->next = p2;p2->last = p3;p1->next = p5;p5->last = p1;break;
	case 15 : invertedOrder(p3,p4);invertedOrder(p5,p6);				head = p4;p3->next = p1;p1->last = p3;p2->next = p6;p6->last = p2;break;
	case 16 : invertedOrder(p1,p2);invertedOrder(p3,p4);invertedOrder(p5,p6);	head = p4;p3->next = p2;p2->last = p3;p1->next = p6;p6->last = p1;break;

	case 17 : 									head = p1;p2->next = p5;p5->last = p2;p6->next = p3;p3->last = p6;break;
	case 18 : invertedOrder(p5,p6);							head = p1;p2->next = p6;p6->last = p2;p5->next = p3;p3->last = p5;break;
	case 19 : invertedOrder(p3,p4);							head = p1;p2->next = p5;p5->last = p2;p6->next = p4;p4->last = p6;break;
	case 20 : invertedOrder(p3,p4);invertedOrder(p5,p6);				head = p1;p2->next = p6;p6->last = p2;p5->next = p4;p4->last = p5;break;

	case 21 : invertedOrder(p1,p2);							head = p2;p1->next = p5;p5->last = p1;p6->next = p3;p3->last = p6;break;
	case 22 : invertedOrder(p1,p2);invertedOrder(p5,p6);				head = p2;p1->next = p6;p6->last = p1;p5->next = p3;p3->last = p5;break;
	case 23 : invertedOrder(p1,p2);invertedOrder(p3,p4);				head = p2;p1->next = p5;p5->last = p1;p6->next = p4;p4->last = p6;break;
	case 24 : invertedOrder(p1,p2);invertedOrder(p3,p4);invertedOrder(p5,p6); 	head = p2;p1->next = p6;p6->last = p1;p5->next = p4;p4->last = p5;break;
	default : printf("type must be more than 1 and less than 25.\n");
    }
    return head;
}

/*
2-opt迭代算法优化path初始解
每次迭代选取两个位点Pi,Pj；如果改变Pi和Pj的节点path指向，得到的路径具有最小的权值，则保留新路径。
重复判断优化，直到路径具有的权值最小
*/
Path * TwoOpt(Path *head,int markNum,float *dist,int *flag_p)
{
    float (*r)[markNum] = (float (*) [markNum])dist;
    int flag;//0表示重复终止，1表示还需要继续重复

    do{
	flag = 0;

	Path *tmp = head;
	Path *end;
	int ftype = 1;//最终确定的重组类型标号
	Path *pi,*py;//遍历Path
	Path *p1,*p2,*p3,*p4,*p5,*p6;//标记重组位置
	Path *fp1,*fp2,*fp3,*fp4,*fp5,*fp6;//最终确定的重组位置
	int x1,x2,x3,x4,x5,x6;//mark编号
	int i,j;
	float finalMinChange = MAXCHANGE;

	while(tmp){if(tmp->next == NULL){end = tmp;}tmp = tmp->next;}//找到end末节点

	pi = head->next->next;
	p1 = head;
	x1 = p1->markIndex;
	p6 = end;
	x6 = p6->markIndex;
	
	extern int UselessMark;
	markNum = markNum - UselessMark;
	for(i = 0;i < markNum - 6;i++)//成三部分，所以p2不能大于markNum - 4,每一个部分的节点数不能小于2
	{
	    p2 = pi;
	    x2 = p2->markIndex;
	    p3 = pi->next;
	    x3 = p3->markIndex;
	    py = pi->next->next;
	    for(j = i + 2;j < markNum - 6;j++)//分成三部分，所以p2不能指向末端,因此j < markNum - 2
	    {
	        int type = 1;//重组类型编号
		p4 = py;
		p5 = py->next;
		x4 = p4->markIndex;
		x5 = p5->markIndex;


		float D_13 = *(*(r + x1) + x3); 
		float D_14 = *(*(r + x1) + x4); 
		float D_15 = *(*(r + x1) + x5); 
		float D_16 = *(*(r + x1) + x6); 
		float D_23 = *(*(r + x2) + x3); 
		float D_24 = *(*(r + x2) + x4); 
		float D_25 = *(*(r + x2) + x5); 
		float D_26 = *(*(r + x2) + x6); 
		float D_35 = *(*(r + x3) + x5); 
		float D_36 = *(*(r + x3) + x6); 
		float D_45 = *(*(r + x4) + x5); 
		float D_46 = *(*(r + x4) + x6); 

		float nochange = D_23 + D_45;// 2 3, 4 5
		float minchange = nochange;//循环第一次时附上初始值

		if(nochange > 0)
		{
		    float change;
		    //20151029改进算法，添加判断条件
		    change = D_24 + D_35;if(change < minchange && change > 0 && D_24 <= D_13 && D_35 <= D_46){minchange = change;type = 2;}//24,35
		    change = D_23 + D_46;if(change < minchange && change > 0 && D_23 <= D_14 && D_46 <= D_35 && D_15 > D_16){minchange = change;type = 3;}//23,46
		    change = D_24 + D_36;if(change < minchange && change > 0 && D_24 <= D_13 && D_36 <= D_45 && D_15 > D_16){minchange = change;type = 4;}//24,36

		    change = D_13 + D_45;if(change < minchange && change > 0 && D_13 <= D_24 && D_45 <= D_36 && D_26 > D_16){minchange = change;type = 5;}//13,45
		    change = D_14 + D_35;if(change < minchange && change > 0 && D_14 <= D_23 && D_35 <= D_46 && D_26 > D_16){minchange = change;type = 6;}//14,35
		    change = D_13 + D_46;if(change < minchange && change > 0 && D_13 <= D_24 && D_46 <= D_35 && D_25 > D_16){minchange = change;type = 7;}//13,46
		    change = D_14 + D_36;if(change < minchange && change > 0 && D_14 <= D_23 && D_36 <= D_45 && D_25 > D_16){minchange = change;type = 8;}//14,36

		    change = D_14 + D_25;if(change < minchange && change > 0 && D_14 <= D_23 && D_25 <= D_16 && D_36 > D_16){minchange = change;type = 9;}//41,25
		    change = D_24 + D_15;if(change < minchange && change > 0 && D_24 <= D_13 && D_15 <= D_26 && D_36 > D_16){minchange = change;type = 10;}//42,15
		    change = D_14 + D_26;if(change < minchange && change > 0 && D_14 <= D_23 && D_26 <= D_15 && D_35 > D_16){minchange = change;type = 11;}//41,26
		    change = D_24 + D_16;if(change < minchange && change > 0 && D_24 <= D_13 && D_16 <= D_25 && D_35 > D_16){minchange = change;type = 12;}//42,16

		    change = D_13 + D_25;if(change < minchange && change > 0 && D_13 <= D_24 && D_25 <= D_16 && D_46 > D_16){minchange = change;type = 13;}//31,25
		    change = D_23 + D_15;if(change < minchange && change > 0 && D_23 <= D_14 && D_15 <= D_26 && D_46 > D_16){minchange = change;type = 14;}//32,15
		    change = D_13 + D_26;if(change < minchange && change > 0 && D_13 <= D_24 && D_26 <= D_15 && D_45 > D_16){minchange = change;type = 15;}//31,26
		    change = D_23 + D_16;if(change < minchange && change > 0 && D_23 <= D_14 && D_16 <= D_25 && D_45 > D_16){minchange = change;type = 16;}//32,16

		    change = D_25 + D_36;if(change < minchange && change > 0 && D_25 <= D_16 && D_36 <= D_45 && D_14 > D_16){minchange = change;type = 17;}//25,63
		    change = D_26 + D_35;if(change < minchange && change > 0 && D_26 <= D_15 && D_35 <= D_46 && D_14 > D_16){minchange = change;type = 18;}//26,53
		    change = D_25 + D_46;if(change < minchange && change > 0 && D_25 <= D_16 && D_46 <= D_35 && D_13 > D_16){minchange = change;type = 19;}//25,64
		    change = D_26 + D_45;if(change < minchange && change > 0 && D_26 <= D_15 && D_45 <= D_36 && D_13 > D_16){minchange = change;type = 20;}//26,54

		    change = D_15 + D_36;if(change < minchange && change > 0 && D_15 <= D_26 && D_36 <= D_45 && D_24 > D_16){minchange = change;type = 21;}//15,63
		    change = D_16 + D_35;if(change < minchange && change > 0 && D_16 <= D_25 && D_35 <= D_46 && D_24 > D_16){minchange = change;type = 22;}//16,53
		    change = D_15 + D_46;if(change < minchange && change > 0 && D_15 <= D_26 && D_46 <= D_35 && D_23 > D_16){minchange = change;type = 23;}//15,64
    		    change = D_16 + D_45;if(change < minchange && change > 0 && D_16 <= D_25 && D_45 <= D_36 && D_23 > D_16){minchange = change;type = 24;}//16,54
		
		    if(minchange < finalMinChange && type > 1)//type大于1，是判断是否需要重组
		    {
		        finalMinChange = minchange;
		        ftype = type;
		        fp1 = p1;
		        fp2 = p2;
		        fp3 = p3;
		        fp4 = p4;
		        fp5 = p5;
		        fp6 = p6;
		    }
		}
		py = py->next;
	    }
	    pi = pi->next;
	}

	if(ftype > 1)//重组
	{
	    head = changePath(fp1,fp2,fp3,fp4,fp5,fp6,ftype);
	    flag = 1;
	    *flag_p = 1;//*flag_p用于main函数，是否继续迭代
	}
    }while(flag);
return head;
}

//根据pi和pj重定位path，pj为NULL时，表示在path头部插入
Path * relocateNode(Path *path,Path *pi,Path *pj)
{
    Path * head = path;
    if(!pj)//在path头部插入
    {
	//连接p(i-1) 和 p(i+1)
	pi->last->next = pi->next;
	if(pi->next){pi->next->last = pi->last;}//pi不是末节点
	
	//pi连接到头部
	head->last = pi;
	pi->next = head;
	pi->last = NULL;
	return pi;
    }else//在Path中部或尾部插入
    {
	//连接p(i-1) 和 p(i+1)
	if(pi->last)//当pi不是头结点
	{
	    pi->last->next = pi->next;
	    if(pi->next)//pi不是末节点，如果是末节点，则无需执行此步
	    {
		pi->next->last = pi->last;
	    }
	}else//当pi是头结点时，则head变成pi->next
	{
	    if(pi->next)
	    {
		head = pi->next;
		head->last = NULL;
	    }
	}

	//连接p(i) 和 p(j+1)
	if(pj->next)
	{
	    pj->next->last = pi;
	    pi->next = pj->next;
	}else
	{
	    pi->next = NULL;
	}
	
	//连接p(i-1) 和 p(i+1)
	//连接p(i) 和 p(j)
	pj->next = pi;
	pi->last = pj;
	return head;
    }
}

/*
我们尝试重定位路径上每个节点到其他所有可能的位置。如果重定位减少了权重，则保留新路径
*/
Path * relocatePath(Path *path,int markNum,float *dist,int *flag_p)
{
    float (*r)[markNum] = (float (*) [markNum])dist;
    float minchange;
    Path *head = path;
    int flag_t;
    do
    {
	flag_t = 0;
        minchange = 0;
        Path *p1,*p2;//遍历Path
        Path *pi,*pj;//标记重组位置,将pi插入pj的后面
        p1 = head;
        while(p1)
        {
	    //当p1节点在path的头部插入时的情况
	    double change;//change = dist(i,head) + dist(i-1,i+1) - dist(i-1,i) - dist(i,i+1)
	    if(p1->last)//当p1不是指向head时
	    {
		if(p1->next)
		{
		    change = *(*(r + p1->markIndex) + head->markIndex);//不能写在一起，否则会造成精度错误
		    change += *(*(r + p1->last->markIndex) + p1->next->markIndex);
		    change -= *(*(r + p1->markIndex) + p1->next->markIndex);
		    change -= *(*(r + p1->markIndex) + p1->last->markIndex);
		    if(*(*(r + p1->markIndex) + head->markIndex) == 0 || *(*(r + p1->last->markIndex) + p1->next->markIndex) == 0){change = 0;}//当重组距离为0时，跳过
		}else
		{
		    change = *(*(r + p1->markIndex) + head->markIndex);
		    change -= *(*(r + p1->markIndex) + p1->last->markIndex);
		    if(*(*(r + p1->markIndex) + head->markIndex) == 0 ){change = 0;}//当重组距离为0时，跳过
		}
	        if((fabs(change) >= 0.000001) && change < minchange)
	        {
		    //printf("minchange=%lf\n",change);
	    	    minchange = change;
		    pi = p1;
		    pj = NULL;
	        }
	    }
            p2 = head;//此初表名mark数必须3
            while(p2)
            {
		if(p1 == p2){p2 = p2->next;continue;}
		float dist1,dist2,dist3,dist4,dist5,dist6;
		//重组前权值
		if(p1->last){dist1 = *(*(r + p1->last->markIndex) + p1->markIndex);}else{dist1 = 0;}//dist1 = dist(i-1,i)
		if(p1->next){dist2 = *(*(r + p1->markIndex) + p1->next->markIndex);}else{dist2 = 0;}//dist2 = dist(i,i+1)
		if(p2->next){dist3 = *(*(r + p2->markIndex) + p2->next->markIndex);}else{dist3 = 0;}//dist3 = dist(j,j+1)

		//重组后权值
		//dist4 = dist(i-1,i+1)
		if(p1->last && p1->next)
		{
		    dist4 = *(*(r + p1->last->markIndex)+ p1->next->markIndex);
		}else
		{
		    dist4 = 0;
		}

		//dist5 = dist(i,j)
		dist5 = *(*(r + p2->markIndex)+ p1->markIndex);
		
		//dist6 = dist(i,j+1)
		if(p2->next)
		{
		    dist6 = *(*(r + p1->markIndex) + p2->next->markIndex);
		}else
		{
		    dist6 = 0;
		}

		if(dist4 == 1000 || dist5 == 1000 || dist6 == 1000){p2 = p2->next;continue;}
//		if((p1->next && p1->next == p2) || (p2->next && p2->next == p1) ){//当p1和p2相邻时，避免精度导致的计算错误
//		float before = dist1 + dist2 + dist3;
//		float after = dist4 + dist5 + dist6;
//		change = after - before;
//		printf("%f,%f,%f\n",after,before,change);
//		if(dist1 == dist4 && dist2 == dist5 && dist3 == dist6){change = 0;}
//		else{change = dist4 + dist5 + dist6 - dist1 - dist2 - dist3;}
//		printf("%f,%f,%f\n%f,%f,%f\n",dist1,dist2,dist3,dist4,dist5,dist6);
//		}else{
//		    if(dist4 == dist2 && dist5 == dist3 && dist6 == dist1){change = 0;}//当p1和p2距离为1时，避免精度导致的计算错误
//		    else if(dist1 == dist4 && dist2 == dist5 && dist3 == dist6){change = 0;}//其他情况，避免精度导致的计算错误
//		    else{change = dist4 + dist5 + dist6 - dist1 - dist2 - dist3;}//change = dist(i-1,i+1) + dist(i,j) + dist(i,j+1) - dist(i-1,i) - dist(i,i+1) - dist(j,j+1)
//		}

		double before = dist1;
		before += dist2;
		before += dist3;
		double after = dist4;
		after += dist5;
		after += dist6;
		change = after - before;
		//printf("%f,%f,%f\n",after,before,change);		
		if(fabs(change) < 0.000001){p2 = p2->next;;continue;}//此处可设置排序的精度

                if(change < minchange)
                {
		    //printf("%f,%f,%f,%f,%f,%f\n",dist1,dist2,dist3,dist4,dist5,dist6);
                    minchange = change;
	  	    //printf("%lf,%lf,%f\n",after,before,change);		
  		    //printf("minchange=%f\n",minchange);
                    pi = p1;
                    pj = p2;
                }
                p2 = p2->next;
            }
            p1 = p1->next;
        }
        if(minchange < 0){	
		//printf("%f,%d,%d,%d,%d\n",minchange,pi->markIndex,pi->next->markIndex,pj->markIndex,pj->next->markIndex);
		head = relocateNode(head,pi,pj);
		*flag_p = 1;
		flag_t = 1; 
	}//保留新路径 *flag_p用于main函数，是否继续迭代
    }while(flag_t);
    return head;
}

//把xls格式转换成SVG格式
void xls2svg(char *title,char *OrderOutputXLS,char *OrderOutputSVG,char *OrderOutputBIN,int markNum){
	
//-------------获取标题，去除路径
	int titleLen = strlen(title);
	int t1 = 0;
	int t2 = 0;
	char finalTitle[titleLen];
	while(t1 < titleLen){
		if(title[t1] == '/' || title[t2] == '\\'){
			strcpy(finalTitle," ");
			t2 = 0;
		}else{
			finalTitle[t2] = title[t1];
			t2++;
		}
		t1++;
	}
	finalTitle[t2] = '\0';

//----------------开始画图
        FILE *infile;
	FILE *binfile;
        FILE *outfile;
        infile = fopen(OrderOutputXLS,"r");
        outfile = fopen(OrderOutputSVG,"w");
        binfile = fopen(OrderOutputBIN,"w");

        int finalMarkNum = markNum - UselessMark;//UselessMark是全局变量
        float lastDist[finalMarkNum];
        float headDist[finalMarkNum];
        char markName[finalMarkNum][30];

	int i = 0;
        int bin = 0;//bin的数量
        int flag = 1;
        char binInfor[10000000] = "";//bin的信息
        while(!feof(infile)){
                char name[30];
                float lastdist;
                float headdist;
                fscanf(infile,"%s\t%f\t%f\n",name,&lastdist,&headdist);
                if(lastdist == 0 && i != 0){
                        i--;
                        if(flag){
                                flag = 0;//下一个mark如果还是同属一个bin，不需要在dist距离上操作，只需要在保存mark的信息
                                bin++;
                                char tempbinNum[5] = "0";
                                char tempbin[10] = "bin";//bin
                                sprintf(tempbinNum,"%d",bin);//将数据格式化的转换成字符串格式
                                strcat(tempbin,tempbinNum);//bin1

                                strcat(binInfor,"\n");
                                strcat(binInfor,tempbin);
                                strcat(binInfor,"\t");
                                strcat(binInfor,markName[i]);
                                strcat(binInfor,",");
                                strcat(binInfor,name);

                                strcpy(markName[i],tempbin);

                        }else{
                                strcat(binInfor,",");
                                strcat(binInfor,name);
                        }
                }else{
                        strcpy(markName[i],name);
                        lastDist[i] = lastdist;
                        headDist[i] = headdist;
                        flag = 1;
                }
                i++;
        }
        fclose(infile);

	double svgHeight = headDist[finalMarkNum - 1] * 5 + 200;
	double chrHeight = headDist[finalMarkNum - 1] * 5 + 20;

        fprintf(outfile,"<?xml version=\"1.0\"?>\n");
        fprintf(outfile,"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n");
        fprintf(outfile,"<svg width=\"500\" height=\"%f\" xmlns=\"http://www.w3.org/2000/svg\">\n",svgHeight);
	fprintf(outfile,"<text x=\"130\" y=\"40\" font-family=\"Verdana\" font-size=\"20\" fill=\"bleak\" >%s</text>",finalTitle);
	fprintf(outfile,"<text x=\"130\" y=\"50\" font-family=\"Verdana\" font-size=\"10\" fill=\"bleak\" >makes num: %d</text>",markNum);
	fprintf(outfile,"<text x=\"130\" y=\"60\" font-family=\"Verdana\" font-size=\"10\" fill=\"bleak\" > bins num: %d</text>",bin);
        fprintf(outfile,"<rect x=\"200\" y=\"90\" width=\"10\" rx=\"4\" ry=\"4\" height=\"%f\" fill=\"red\" stroke=\"red\"/>\n",chrHeight);//染色体长条

        flag = 100;//上一个mark的位置
        float last_loc = 0;
	i = 0;
        while(i < finalMarkNum){
                float headdist = headDist[i];

                last_loc = headdist;
                float text_y = last_loc;//文字y坐标
                float heng_y = text_y * 5 + 100;
                if(heng_y >= flag){//该位点没有被占据
                	text_y = heng_y + 4;
                        flag = text_y + 10;
                }else{
                        text_y = flag;
                        flag = text_y + 10;
                }
                float zhe_y = text_y - 4;

                fprintf(outfile,"\n<text x=\"157\" y=\"%f\" font-family=\"Verdana\" font-size=\"10\" text-anchor = \"end\" fill=\"bleak\" >%5.2f</text>\n",text_y,last_loc);
                fprintf(outfile,"<polyline  points=\"158,%f 200,%f 210,%f 260,%f\"  stroke=\"blue\" stroke-width=\"1\" fill=\"none\"/>\n",zhe_y,heng_y,heng_y,zhe_y);
                fprintf(outfile,"<text x=\"260\" y=\"%f\" font-family=\"Verdana\" font-size=\"10\" fill=\"bleak\" >%s</text>\n\n",text_y,markName[i]);

		i++;
        }
	fprintf(binfile,"%s\n",binInfor);
        fprintf(outfile,"</svg>\n");
        fclose(outfile);
        fclose(binfile);

}

//测试用
void printDist(float *dist,int markNum){
	float (*r)[markNum] = (float (*) [markNum])dist;
	int i = 0;
	int j = 0;
	for(i=0;i<markNum;i++){
		for(j=0;j<markNum;j++){
			printf("%f\t",*(*(r+i)+j));
		}
		printf("\n");
	}
}

int step3_order(char *infile,int function)
{

    fprintf(stdout,"********************\n");
    fprintf(stdout,"Order\n");
    fprintf(stdout,"********************\n");
    fprintf(stdout,"Parameters: u2map order -f3 %s -f %d\n",infile,function);

    if(infile == NULL){
	fprintf(stderr,"Usage:u2map order -f3 MLE.xls\n");
	fprintf(stderr,"\n\t-f3\t[str] input file for ordering.\n");
	fprintf(stderr,"\t-f\t[int] mapping function.1 for kosambi, 2 for haldane. [defalut 1]\n\n");
	exit(EXIT_FAILURE);

    }

    //读取mark数量
    FILE *fp;
    int markNum = 0;//MARK的数量
    if((fp = fopen(infile,"r")) != NULL){ 
        fscanf(fp,"markNum=%d",&markNum);
        fclose(fp);
    }else{
	fprintf(stderr,"ERROR: can not open %s file\n",infile);
	exit(EXIT_FAILURE);
    }
    if(markNum < 6){
	fprintf(stderr,"Warn: %s mark numbes less than 6,skipping...\n",infile);
	return 0;
    }
    //输出xls&svg文件
    int dirLen = strlen(infile) + 20;
    char OrderOutputXLS[dirLen];
    char OrderOutputBIN[dirLen];
    char OrderOutputSVG[dirLen];
    strcpy(OrderOutputXLS,infile);
    strcpy(OrderOutputBIN,infile);
    strcpy(OrderOutputSVG,infile);
    strcat(OrderOutputXLS,".order.xls");
    strcat(OrderOutputBIN,".bin_infor.list");
    strcat(OrderOutputSVG,".order.svg");

    //#######step1########
    //读取markName和r[][]重组率矩阵，并用Kosambi作图函数更好地表示重组率
    char markName[markNum][30];
    float r[markNum][markNum];
    float *dist = (float *)r;
    readGroupTable(infile,dist,markName,function);
    BinNode *useless = NULL;
    useless = mergeTable(dist,markNum);
    if(markNum - UselessMark < 6){fprintf(stderr,"Warn: %s useful mark numbes less than 6,skipping...\n",infile);return 0;}
    //printf("UselessMark = %d\n",UselessMark);
    //printDist(dist,markNum);
    //构建最小生成树
    MARKNode *root = (MARKNode *)malloc(sizeof(MARKNode));
    MST(dist,markNum,root);
    //获取最小生成树中最长的path，以及生成树以外的其他节点串
    modifyLevel(root);//更新最小生成树的level值
    //printMST(root,markName);//输出最小生成树
    getLongestPathRoot(root);//找到最长路径的root节点
    //printMST(pathRoot,markName);
    Path *orderPath = selectLongestPath(pathRoot);//orderPath最长路径的头结点,pathRoot是全局变量，在函数getLongestPathRoot中赋值
    //printPath(orderPath,markName,dist,markNum);
    Path *noOrderPath = getNoOrder(orderPath,markNum,dist);//noOrderPath尚未在orderPath中的mark节点
    //printf("\n*******MST*******\n");
    //printf("\n*******Not in MST*******\n");
    //printPath(noOrderPath,markName);    
    //#######step2########
    //将noOrderPath的节点整合到orderPath，如果已经是一条完整的path，则跳过此步骤
    readGroupTable(infile,dist,markName,function);//更新r[][]矩阵，因为在建树过程中有修改其中的值

    orderPath = mergePath(orderPath,noOrderPath,dist,markNum);
    //printf("\n*******After Merge*******\n");
    //printPath(orderPath,markName,dist,markNum);

    //我们迭代地应用2-opt 启发法和重定位启发法直到都不能进一步改良路径权重。所得TSP 路径代表了我们的最终解.
    //printf("\n*******TwoOpt*******\n");
    int flag = 1;
    int *flag_p = &flag;
    orderPath = TwoOpt(orderPath,markNum,dist,flag_p);
    flag = 1;//先进行一次2-opt 改良，不管有无改进，都需要进行重定位改良
    int i = 0;
    while(flag)
    {
        //printf("markNum = %d,i = %d\n",markNum,i);
	flag = 0;
	if(i%2){orderPath = TwoOpt(orderPath,markNum,dist,flag_p);}//flag_p会更新
	else{orderPath = relocatePath(orderPath,markNum,dist,flag_p);}
	i++;
	if(i>100){fprintf(stderr,"Warn: repeat to run TwoOpt and relocatePath more than 100 times.\n");break;}//一般都是重复几次，当迭代次数太多时，代表进入死循环。人工跳出
    }
    printPath(orderPath,markName,dist,markNum,useless,OrderOutputXLS);
    xls2svg(infile,OrderOutputXLS,OrderOutputSVG,OrderOutputBIN,markNum);
    fprintf(stdout,"Order finish!\n\n");
    return 0;
}

