#ifndef _ORDER_H_
#define _ORDER_H_

#define MAXRLE 64425
#define MAXCHANGE 2000
#define UNLESS 2000

typedef struct markNode{
    int markIndex;//markName 的下标
    int Level;//最长孩子节点的层数，最下层值为0
    struct markNode *parent;
    struct markNode *leftChild;
    struct markNode *rightSibling;
}MARKNode;//相当于图的顶点

typedef struct treeList{
    int markIndex;//markName 的下标
    struct markNode *mstNode;
    struct treeList *next;//下一个节点
}TreeList;//用于存储已建树的MARK 下标

typedef struct PathNode{
    int markIndex;
    struct PathNode *next;
    struct PathNode *last;
}Path;//用于存储最长MST的节点和其他分支节点

typedef struct binNode{
        int keepMark;//保留的marks下标
        struct PathNode *otherMarks;//与之距离为0的marks下标
        struct binNode *nextBin;
}BinNode;//在构建MST树之前，将重组率为0的marks整合成一个bin

void MST(float *dist,int markNUM,MARKNode *root);
float MLE2MapingDistance(float r,int function);
void readGroupTable(char file[],float *r,char (*name)[30],int function);
BinNode * mergeTable(float *dist,int markNum);
int modifyLevel(MARKNode *node);
int printMST(MARKNode *node,char (*name)[30]);
MARKNode * maxChild(MARKNode *node);
void getLongestPathRoot(MARKNode *root);
Path * getNoOrder(Path *orderPath,int markNum,float *dist);
Path * selectLongestPath(MARKNode *root);
BinNode * isBin(int index,BinNode *head,float distance,char (*name)[30],char *OrderOutputXLS);
void printPath(Path *path,char (*name)[30],float *dist,int markNum,BinNode *useless,char *OrderOutputXLS);
Path * mergePath(Path *orderHead,Path *noOrderHead,float *dist,int markNum);
void invertedOrder(Path *p1,Path *p2);
Path * changePath(Path *p1,Path *p2,Path *p3,Path *p4,Path *p5,Path *p6,int type);
Path * TwoOpt(Path *head,int markNum,float *dist,int *flag_p);
Path * relocateNode(Path *path,Path *pi,Path *pj);
Path * relocatePath(Path *path,int markNum,float *dist,int *flag_p);
void xls2svg(char *title,char *OrderOutputXLS,char *OrderOutputSVG,char *OrderOutputBIN,int markNum);
void printDist(float *dist,int markNum);
int step3_order(char *infile,int function);

#endif
