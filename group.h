#ifndef _GROUP_H_
#define _GROUP_H_

#define MAX 999999
typedef struct SimilarNode{
    struct ClusterNode *clu;//指向相似的簇
    struct SimilarNode *next;//下一个
}Similar;

typedef struct MarkNode{
    int markno;//基因标记标号
    struct MarkNode *nextMark;//下一个标记
}Mark;

typedef struct ClusterNode{
    int clusterNo;//簇编号
    int num;//该簇所包含的MARK数量
    struct MarkNode *mark;//该簇所包含的MARK链表
    struct SimilarNode *similar;//距离最近的簇指针列表
    float similarDistance;//最近的簇的距离
    int similarNum;//距离最近的簇的个数
    struct ClusterNode *lastCluster;//上一个簇
    struct ClusterNode *nextCluster;//下一个簇
}Cluster;

typedef struct MinDist{
//第n次聚类的簇列表(20150611修改 新增last属性)
    struct MinDist *last;
    struct MinDist *next;
    struct SimilarNode *sim;
}MinDistList;

//聚类过程距离记录链表，用于输出聚类过程的信息，供用户参考筛选合适的分群
typedef struct DistList{
    float dist;//当前聚类最小距离
    int cluNum;//剩余簇的数量
    struct DistList *next;//下一次聚类描述
}ClusterDist;

void InitCluster(Cluster *clus,int MARK);
float distMin(Cluster *cluster1,Cluster *cluster2,float *r,int MARK);
float distMax(Cluster *cluster1,Cluster *cluster2,float *r,int MARK);
float distAver(Cluster *cluster1,Cluster *cluster2,float *r,int MARK);
void freeSimilar(Similar * si);
int updateCluster(Cluster *clus,float *r,int MARK,int method);
void freeDistList(MinDistList *minDist);
float SelectCluster(MinDistList *minDist,Cluster *clus,int method);
int isRelated(Similar *p1,Similar *p2);
void mergeSimilar(Similar *p1,Similar *p2);
void mergeClusterList(MinDistList *minDist);
Cluster *mergeCluster(MinDistList *minDist,Cluster *clus);
void printCluster(Cluster * clus,char filename1[1000],char filename2[1000],char (*name)[30],float distSubtraction);
void readTable(char file[],float *r,char (*name)[30]);
void readLODTable(char file[],float *r);
void printToFile(char *inputfile,char *outputDir,float *dist,char (*name)[30], int MARK,float index,int order,int function,int method);
void GroupErrorPrint();
int step2_group(char *file,char *output,int CHR,int TOP_Num,int order,int function,int method,float lod1,float lod2,char *lodfile,float step);

#endif
