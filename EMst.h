#pragma once

#include<iostream>
#include<string>
#include<stdio.h>
#include<algorithm>
#include<cmath>
#include <fstream>
#include <vector>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/common_headers.h>
#include <pcl/kdtree/kdtree_flann.h>

#define path "G:/Code/new/20230801/k-mst/k-mst/lookprocess"
//"E:/PointCloud/Code/k-mst/k-mst/lookprocess/" 1014

typedef int VertexType;
typedef float EdgeType;
#define MAXVEX 5000

#define Inf 255800.000
using namespace std;

//the structure of data in MST
typedef struct EdgeNode {
	int adjvex;
	EdgeType weight;
	struct EdgeNode *next;
}EdgeNode;

typedef struct VertexNode {
	VertexType data;
	int degree;
	EdgeNode *firstedge;
	bool visit;
	bool is_cor;
}AdjList;

typedef struct {
	AdjList *adjList;
	int numVertsxes, numEdges;
}GraphAdjList;

struct EndNode {
	int data;
	int degree;
	vector<int> lineorder;
};

typedef struct degreeMore2Node
{
	int data;
	int degree;
	vector<EndNode> endnode;
	int maxlength;

}crticalNode;

struct save {
	int depthnum;
	vector<int> dirdata;
};

struct closedge {
	float lowcost; //closedge[i].lowcost表示当前与第i个顶点相连的的权值最小的边
	int vex; //closedeg[i].vex表示与第i个顶点相连的顶点的位置 
};
// the structure of point xyzp
struct MyPointType
{
	PCL_ADD_POINT4D;  //该点类型有4个元素      

	/*尝试新增一个自定义*/
	float intensity;
	float gvalue;
	float p0;
	float p1;
	float ptype;


	EIGEN_MAKE_ALIGNED_OPERATOR_NEW   //确保new操作符对齐操作 

}EIGEN_ALIGN16;   //强制SSE 对齐

POINT_CLOUD_REGISTER_POINT_STRUCT(MyPointType,    //注册点类型宏
(float, x, x)
(float, y, y)
(float, z, z)
(float, intensity, intensity)
(float, gvalue, gvalue)
(float, p0, p0)
(float, p1, p1)
(float, ptype, ptype)
)

struct RawPointType    //定义点类型结构
{
	PCL_ADD_POINT4D;  //该点类型有4个元素      

	/*尝试新增一个自定义*/
	float intensity;
	float gvalue;
	float gx;
	float gy;
	float gz;
	float Tlambda0;
	float Tlambda1;
	float Tlambda2;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW   //确保new操作符对齐操作 

}EIGEN_ALIGN16;   //强制SSE 对齐

POINT_CLOUD_REGISTER_POINT_STRUCT(RawPointType,    //注册点类型宏
(float, x, x)
(float, y, y)
(float, z, z)
(float, intensity, intensity)
(float, gvalue, gvalue)
(float, gx, gx)
(float, gy, gy)
(float, gz, gz)
(float, Tlambda0, Tlambda0)
(float, Tlambda1, Tlambda1)
(float, Tlambda2, Tlambda2)
)

typedef pcl::PointCloud<pcl::PointXYZ> PointXYZ;
typedef pcl::PointCloud<pcl::PointXYZ>::Ptr p_PointXYZ;

typedef pcl::PointCloud<MyPointType> Mypoint;
typedef pcl::PointCloud<MyPointType>::Ptr  p_Mypoint;

typedef pcl::PointCloud<RawPointType> rawPoint;
typedef pcl::PointCloud<RawPointType>::Ptr  p_rawPoint;

void pruning4_plus(GraphAdjList* tree1, int index, int* depth, vector<int> &a);
static EdgeNode * make_node(const int pos, const float distance);
static void initial_graph(GraphAdjList * graph, GraphAdjList * kruskal_tree);

template<typename T, typename... U>
void logger(T t, U... ts);


///// 声明粒子滤波函数
void particleFilter();

///// 声明获取权值函数
	//1018***********//float getWeightValue(int i, int j); // 添加 getWeightValue 函数声明
float getWeightValue(int i, int j);


class EMst
{
public:
	GraphAdjList Graph;
	//Mypoint点云数据类型为5维<I,g,p0,p1,ptype>,该数据
	Mypoint rawpoint;
	//
	Mypoint knnPoint;

	//rawPoint点云类型为8维数据<I,g,gx,gy,gz,lamb1,lamb2,lamb3>
	rawPoint raw;


public:
	EMst(string datapath, bool flag);
	void CreateALGraph(p_Mypoint P, int K)//创建无向图结构 1014
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::copyPointCloud(*P, *cloud);
		pcl::KdTreeFLANN<pcl::PointXYZ> searchtree;
		searchtree.setInputCloud(cloud);
		pcl::PointXYZ searchpoint;
		std::vector<int> Neighbors;
		std::vector<float> SquaredDistance;
		int i, j, k;
		EdgeNode* e;
		cout << "knn点数：" << this->Graph.numVertsxes << endl;


		for (i = 0; i < this->Graph.numVertsxes; i++) {
			//ofstream fp("po.txt", ios::app);
			//fp << P->points[i].x << " " << P->points[i].y << " " << P->points[i].z << " " << P->points[i].p0 << " " << P->points[i].p1 << " " << endl;
			searchpoint.x = P->points[i].x;
			searchpoint.y = P->points[i].y;
			searchpoint.z = P->points[i].z;
			searchtree.nearestKSearch(searchpoint, K, Neighbors, SquaredDistance);
			float pi = P->points[i].p0;  //当前点的概率
			float pj = 0.0;
			float wij = 0.0;
			// 权值函数（粒子滤波）
			for (k = 1; k < Neighbors.size(); k++) {
				j = Neighbors[k]; //全局索引号，而k是局部的近邻索引
				pj = P->points[j].p0; //第k个邻域点的概率
				wij = 0; 	//			
				//权值函数 float getWeigtValue(i,j) 引用main.cpp中计算权值的函数，得到wij
				//1018*****************************//float getWeigtValue(i,j)
				e = (EdgeNode*)malloc(sizeof(EdgeNode));
				e->adjvex = j;
				e->weight = wij;//
				e->next = this->Graph.adjList[i].firstedge;
				this->Graph.adjList[i].firstedge = e;

				e = (EdgeNode*)malloc(sizeof(EdgeNode));
				e->adjvex = i;
				e->weight = wij;
				e->next = this->Graph.adjList[j].firstedge;
				this->Graph.adjList[j].firstedge = e;

			}
		
		}
		//for (k = 1; k < Neighbors.size(); k++) {
		//	j = Neighbors[k]; // 全局索引号，而 k 是局部的近邻索引
		//	pj = P->points[j].p0; // 第 k 个邻域点的概率
		//	wij = getWeightValue(i, j); // 使用 getWeightValue 函数计算权值
		//}

		/*for (i = 0; i < this->Graph.numVertsxes; i++) {
			//ofstream fp("po.txt", ios::app);
			//fp << P->points[i].x << " " << P->points[i].y << " " << P->points[i].z << " " << P->points[i].p0 << " " << P->points[i].p1 << " " << endl;
			searchpoint.x = P->points[i].x;
			searchpoint.y = P->points[i].y;
			searchpoint.z = P->points[i].z;
			searchtree.nearestKSearch(searchpoint, K, Neighbors, SquaredDistance);
			float pi = P->points[i].p0;  //当前点的概率
			float pj = 0.0;
			float pij = 0.0;
			float wij = 0.0;
			// 权值函数（粒子滤波）


			for (k = 1; k < Neighbors.size(); k++) {
				if (sqrt(SquaredDistance[k])) {
					j = Neighbors[k];
					pj = P->points[j].p0;
					if (pi == pj) pij = exp(sqrt(SquaredDistance[k])* pi);
					else pij = exp(sqrt(SquaredDistance[k])*(pi * (log(pi) - 1) + pj * (1 - log(pj))) / (pi - pj));
					if (pij >= 1) {
						wij = pi * pj;
					}
				
					else {
						float pijno = exp(sqrt(SquaredDistance[k])*(pi * (log(pi) - 1) + (1-pj) * (1 - log(1-pj))) / (pi - (1-pj)));
						float pinoj = exp(sqrt(SquaredDistance[k])*((1-pi) * (log(1-pi) - 1) + pj * (1 - log(pj))) / ((1-pi) - pj));
						float pinojno = exp(sqrt(SquaredDistance[k])*((1 - pi) * (log(1 - pi) - 1) + (1 - pj) * (1 - log(1-pj))) / ((1 - pi) - (1-pj)));
						wij = pi * pj * -log(pij / (1 - pij)) 
							+  pi * (1 - pj)* -log(pijno / (1 - pijno)) 
							+ (1-pi) * pj * -log(pinoj / (1 - pinoj)) 
							+ (1-pi) * (1 - pj)* -log(pinojno / (1 - pinojno));
					}
					// wij =  -log(pij);

					if (wij < 0) {
						e = (EdgeNode *)malloc(sizeof(EdgeNode));
						e->adjvex = j;
						e->weight = wij;
						e->next = this->Graph.adjList[i].firstedge;
						this->Graph.adjList[i].firstedge = e;

						e = (EdgeNode *)malloc(sizeof(EdgeNode));
						e->adjvex = i;
						e->weight = wij;
						e->next = this->Graph.adjList[j].firstedge;
						this->Graph.adjList[j].firstedge = e;
						
					}
					else {
						e = (EdgeNode *)malloc(sizeof(EdgeNode));
						e->adjvex = j;
						e->weight = 1;
						e->next = this->Graph.adjList[i].firstedge;
						this->Graph.adjList[i].firstedge = e;

						e = (EdgeNode *)malloc(sizeof(EdgeNode));
						e->adjvex = i;
						e->weight = 1;
						e->next = this->Graph.adjList[j].firstedge;
						this->Graph.adjList[j].firstedge = e;
					}

				}
			}
			
		}*/
	}




	void MiniTree_Prim(GraphAdjList G, int v, GraphAdjList * result)
	{
		int i, j;
		bool *visit = new bool[G.numVertsxes];
		closedge *closedges = new closedge[G.numVertsxes];
		initial_graph(&G, result);
		//初始化closedge
		for (i = 0; i < G.numVertsxes; i++)
		{
			visit[i] = false;
			closedges[i].vex = v;
			closedges[i].lowcost = Inf;
		}
		visit[v] = true;
		EdgeNode *p, *tmp;
		p = G.adjList[v].firstedge;
		while (p)
		{
			closedges[p->adjvex].lowcost = p->weight;
			p = p->next;
		}
		//找出closedge中最小的边，输出并进行更新
		for (j = 0; j < G.numVertsxes; j++)
		{
			double min = Inf;
			int t = Inf;
			for (i = 0; i < G.numVertsxes; i++)
			{
				if (closedges[i].lowcost < min && visit[i] == false)
				{
					min = closedges[i].lowcost;
					t = i;
				}
			}
			if (t != Inf) {
				//printf("%d,%d,%f\n", closedges[t].vex, t,closedges[t].lowcost);
				int star = closedges[t].vex;
				int to = t;
				visit[t] = true;
				if (result->adjList[t].firstedge == NULL) {
					result->adjList[t].firstedge = make_node(closedges[t].vex, closedges[t].lowcost);
					if (result->adjList[closedges[t].vex].firstedge == NULL) {
						result->adjList[closedges[t].vex].firstedge = make_node(t, closedges[t].lowcost);
					}
					else {
						tmp = result->adjList[closedges[t].vex].firstedge;
						while (tmp->next != NULL)
							tmp = tmp->next;
						tmp->next = make_node(t, closedges[t].lowcost);
					}
				}
				else {
					tmp = result->adjList[t].firstedge;
					while (tmp->next != NULL)
						tmp = tmp->next;
					tmp->next = make_node(closedges[t].vex, closedges[t].lowcost);
					if (result->adjList[closedges[t].vex].firstedge == NULL) {
						result->adjList[closedges[t].vex].firstedge = make_node(t, closedges[t].lowcost);
					}
					else {
						tmp = result->adjList[closedges[t].vex].firstedge;
						while (tmp->next != NULL)
							tmp = tmp->next;
						tmp->next = make_node(t, closedges[t].lowcost);
					}
				}
				//更新 
				p = G.adjList[t].firstedge;
				while (p)
				{
					if (closedges[p->adjvex].lowcost > p->weight)
					{
						closedges[p->adjvex].lowcost = p->weight;
						closedges[p->adjvex].vex = t;
					}
					p = p->next;
				}

			}
			else {

				logger("The graph whose initial point is ", v, " has been over!");
				//cout << "The graph whose initial point is " << v << " has been over!" << endl;
				break;
			}
		}
		delete[] closedges;
		delete[] visit;
	}
	void findline(GraphAdjList * tree1, int index, vector<EndNode> &line) {
		save *b = new save[tree1->adjList[index].degree];
		if (tree1->adjList[index].degree >= 3)
		{
			int i = 0;
			EdgeNode* p;
			p = tree1->adjList[index].firstedge;
			while (p)
			{
				int depth = 0;
				b[i].dirdata.push_back(index);
				tree1->adjList[index].visit = true;
				pruning4_plus(tree1, p->adjvex, &depth, b[i].dirdata);
				b[i].depthnum = depth;
				++i;
				p = p->next;
			}

		}
		for (int de = 0; de < tree1->adjList[index].degree; de++)
		{
			EndNode temp;
			if (b[de].dirdata.size())
			{
				temp.lineorder = b[de].dirdata;
				//尾结点：度为1 或者 度为大于等于3的结点，因此需要记录下来；
				int end = b[de].dirdata.back();
				temp.degree = tree1->adjList[end].degree;
				temp.data = end;
				line.push_back(temp);
			}

		}

	}
	~EMst()
	{}
};

//
//typedef struct treeNode
//{
//	bool isVisited;
//	int position;
//	treeNode *parent;
//	vector<treeNode *> childList;
//	EdgeNode *vex;
//	vector<EdgeNode *> *seedList;
//}treeNode;

float ProLine = 0.0;
float threshold_p = 0.0;

//flag表示是否取消计算p0，如果为真则不计算
EMst::EMst(string dataPath, bool flag)
{
	//p_Mypoint inputrawdata;
	if (flag){
		pcl::io::loadPCDFile(dataPath, this->rawpoint);
		cout << "点云数据输入完成." << endl;
	}
	else
	{
		pcl::io::loadPCDFile(dataPath, this->raw);
		logger("点云数据输入完成.");
		//cout << "点云数据输入完成." << endl;
		int NumPoint = this->raw.size();

		for (int i = 0; i < NumPoint; i++) {
		/*	float normagvalue = 0.0f;
			if (this->raw.points[i].gvalue > 0.01){
				normagvalue = 1.0f;
			}
			else {
				normagvalue = this->raw.points[i].gvalue * 100;
			} */
			float p0 = 0.5f *(this->raw.points[i].intensity + this->raw.points[i].gvalue);
			//float p0 = 0.5f*(this->raw.points[i].intensity + (1 - exp(-this->raw.points[i].gvalue / 0.2)));
			//logger(p0);
			float p1 = 1 - p0;
			if (p0*p1 < 0)logger("错误发生在第", i, "行");
			else {
				MyPointType p;
				p.x = this->raw.points[i].x;
				p.y = this->raw.points[i].y;
				p.z = this->raw.points[i].z;
				p.p0 = p0;
				p.p1 = p1;
				p.ptype = p0 > 0.5 ? 0 : 1;
				this->rawpoint.points.push_back(p);
			}
		}
	}
	int NumPoint = this->rawpoint.size();

	//define the points in KNN graph

	for (int i = 0; i < NumPoint; i++) {
		ProLine = this->rawpoint.points[i].p0;
		if (ProLine > threshold_p) this->knnPoint.push_back(this->rawpoint.points[i]);
	}
	int NumknnPoint = this->knnPoint.size();
	GraphAdjList *G = &this->Graph;
	G->adjList = new AdjList[NumknnPoint];
	G->numVertsxes = NumknnPoint;
	for (int i = 0; i < NumknnPoint; i++) {
		G->adjList[i].data = i;
		G->adjList[i].degree = 0;
		G->adjList[i].firstedge = NULL;
	}


}


static EdgeNode * make_node(const int pos, const float distance)
{
	EdgeNode * new_node = (EdgeNode *)malloc(sizeof(EdgeNode));
	if (new_node == NULL)
		exit(1);

	new_node->next = NULL;
	new_node->weight = distance;
	new_node->adjvex = pos;

	return new_node;
}

static void initial_graph(GraphAdjList * graph, GraphAdjList * kruskal_tree)
{
	int i;
	kruskal_tree->numVertsxes = graph->numVertsxes;
	kruskal_tree->numEdges = graph->numEdges;

	for (i = 0; i < graph->numVertsxes; i++)
	{
		kruskal_tree->adjList[i].data = graph->adjList[i].data;
		kruskal_tree->adjList[i].firstedge = NULL;
	}
}

void pruning4_plus(GraphAdjList* tree1, int index, int* depth, vector<int> &a)
{
	EdgeNode *p = tree1->adjList[index].firstedge;
	if (tree1->adjList[index].degree == 1)
	{
		tree1->adjList[index].visit = true;
		a.push_back(index);
		(*depth)++;
	}

	if (tree1->adjList[index].degree == 2)
	{
		if (!tree1->adjList[index].visit)
		{
			tree1->adjList[index].visit = true;
			a.push_back(index);
			if (tree1->adjList[p->adjvex].visit)
			{
				if (p->next) {
					p = p->next;
				}
				if (!tree1->adjList[p->adjvex].visit && tree1->adjList[p->adjvex].degree <= 2)
				{
					pruning4_plus(tree1, p->adjvex, depth, a);
				}
				else
				{
					if (tree1->adjList[p->adjvex].degree > 2) {
						a.push_back(p->adjvex);
					}
					else
					{
						p = tree1->adjList[index].firstedge;
						a.push_back(p->adjvex);
					}

				}
			}
			else
			{
				(*depth)++;
				pruning4_plus(tree1, p->adjvex, depth, a);
			}
		}


	}

	if (tree1->adjList[index].degree >= 3)
	{
		tree1->adjList[index].visit = true;
		a.push_back(index);
		(*depth)++;
	}

}
