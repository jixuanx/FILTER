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
	float lowcost; //closedge[i].lowcost��ʾ��ǰ���i�����������ĵ�Ȩֵ��С�ı�
	int vex; //closedeg[i].vex��ʾ���i�����������Ķ����λ�� 
};
// the structure of point xyzp
struct MyPointType
{
	PCL_ADD_POINT4D;  //�õ�������4��Ԫ��      

	/*��������һ���Զ���*/
	float intensity;
	float gvalue;
	float p0;
	float p1;
	float ptype;


	EIGEN_MAKE_ALIGNED_OPERATOR_NEW   //ȷ��new������������� 

}EIGEN_ALIGN16;   //ǿ��SSE ����

POINT_CLOUD_REGISTER_POINT_STRUCT(MyPointType,    //ע������ͺ�
(float, x, x)
(float, y, y)
(float, z, z)
(float, intensity, intensity)
(float, gvalue, gvalue)
(float, p0, p0)
(float, p1, p1)
(float, ptype, ptype)
)

struct RawPointType    //��������ͽṹ
{
	PCL_ADD_POINT4D;  //�õ�������4��Ԫ��      

	/*��������һ���Զ���*/
	float intensity;
	float gvalue;
	float gx;
	float gy;
	float gz;
	float Tlambda0;
	float Tlambda1;
	float Tlambda2;

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW   //ȷ��new������������� 

}EIGEN_ALIGN16;   //ǿ��SSE ����

POINT_CLOUD_REGISTER_POINT_STRUCT(RawPointType,    //ע������ͺ�
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


///// ���������˲�����
void particleFilter();

///// ������ȡȨֵ����
	//1018***********//float getWeightValue(int i, int j); // ��� getWeightValue ��������
float getWeightValue(int i, int j);


class EMst
{
public:
	GraphAdjList Graph;
	//Mypoint������������Ϊ5ά<I,g,p0,p1,ptype>,������
	Mypoint rawpoint;
	//
	Mypoint knnPoint;

	//rawPoint��������Ϊ8ά����<I,g,gx,gy,gz,lamb1,lamb2,lamb3>
	rawPoint raw;


public:
	EMst(string datapath, bool flag);
	void CreateALGraph(p_Mypoint P, int K)//��������ͼ�ṹ 1014
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
		cout << "knn������" << this->Graph.numVertsxes << endl;


		for (i = 0; i < this->Graph.numVertsxes; i++) {
			//ofstream fp("po.txt", ios::app);
			//fp << P->points[i].x << " " << P->points[i].y << " " << P->points[i].z << " " << P->points[i].p0 << " " << P->points[i].p1 << " " << endl;
			searchpoint.x = P->points[i].x;
			searchpoint.y = P->points[i].y;
			searchpoint.z = P->points[i].z;
			searchtree.nearestKSearch(searchpoint, K, Neighbors, SquaredDistance);
			float pi = P->points[i].p0;  //��ǰ��ĸ���
			float pj = 0.0;
			float wij = 0.0;
			// Ȩֵ�����������˲���
			for (k = 1; k < Neighbors.size(); k++) {
				j = Neighbors[k]; //ȫ�������ţ���k�Ǿֲ��Ľ�������
				pj = P->points[j].p0; //��k�������ĸ���
				wij = 0; 	//			
				//Ȩֵ���� float getWeigtValue(i,j) ����main.cpp�м���Ȩֵ�ĺ������õ�wij
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
		//	j = Neighbors[k]; // ȫ�������ţ��� k �Ǿֲ��Ľ�������
		//	pj = P->points[j].p0; // �� k �������ĸ���
		//	wij = getWeightValue(i, j); // ʹ�� getWeightValue ��������Ȩֵ
		//}

		/*for (i = 0; i < this->Graph.numVertsxes; i++) {
			//ofstream fp("po.txt", ios::app);
			//fp << P->points[i].x << " " << P->points[i].y << " " << P->points[i].z << " " << P->points[i].p0 << " " << P->points[i].p1 << " " << endl;
			searchpoint.x = P->points[i].x;
			searchpoint.y = P->points[i].y;
			searchpoint.z = P->points[i].z;
			searchtree.nearestKSearch(searchpoint, K, Neighbors, SquaredDistance);
			float pi = P->points[i].p0;  //��ǰ��ĸ���
			float pj = 0.0;
			float pij = 0.0;
			float wij = 0.0;
			// Ȩֵ�����������˲���


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
		//��ʼ��closedge
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
		//�ҳ�closedge����С�ıߣ���������и���
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
				//���� 
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
				//β��㣺��Ϊ1 ���� ��Ϊ���ڵ���3�Ľ�㣬�����Ҫ��¼������
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

//flag��ʾ�Ƿ�ȡ������p0�����Ϊ���򲻼���
EMst::EMst(string dataPath, bool flag)
{
	//p_Mypoint inputrawdata;
	if (flag){
		pcl::io::loadPCDFile(dataPath, this->rawpoint);
		cout << "���������������." << endl;
	}
	else
	{
		pcl::io::loadPCDFile(dataPath, this->raw);
		logger("���������������.");
		//cout << "���������������." << endl;
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
			if (p0*p1 < 0)logger("�������ڵ�", i, "��");
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
