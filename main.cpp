#include <iostream>
#include<stdio.h>
#include <math.h> 
#include<algorithm>
#include<iostream>
#include<string>
#include <fstream> 
#include <vector> 
#include "helper.h" 
#include "simplify.h"
#include <pcl/point_types.h>
#include <pcl/point_cloud.h> 
#include "EMst.h"

void txtToPCD(string txtname);
int order_element_in_vector(const vector<int> &v, int element);

void mains(EMst &g, int Orderitreator);
//#define path "E:/PointCloud/Code/k-mst/k-mst/result/2_house.pcd.pl"
int main(int argc, char* argv[])
{
	//这里进行循环，计算概率
	//cleanFile(path);

	/*step1: 读取pcd and score*/
	//string pathname = "1_svm_result.txt";
	//txtToPCD(pathname);

	/*step2: 构图，权值为wij*/
	//读取文件 //flag表示是否取消计算p0，如果为真则不计算,如果为假则计算
	//"E:/PointCloud/Data_exp/boundary/fac1-BeforeThinning.pcd"
	//E:/PointCloud/Data_exp/limitation/ALS_BP.pcd

	EMst g("G:/Code/new/20230801/k-mst/k-mst/lookprocess/3c-02.pcd", false);

	for (int i = 0; i < 5; i++)
	{
		cout << "第" << i << "次迭代" << endl;
		mains(g, i);
		clear(deprecatedList);
	}
	//cleanFile(path);
	//mains(g);
}




void mains(EMst &g, int Orderitreator) {

	//PCL点云类型
	pcl::PointCloud<MyPointType>::Ptr m_knnpoint(new pcl::PointCloud<MyPointType>);
	//深度拷贝

	m_knnpoint = (g.knnPoint).makeShared();
	//构图
	g.CreateALGraph(m_knnpoint, 20);

	/////////////////// 创建 EMst 对象
	//////////////EMst g;
	//////////////// 调用粒子滤波算法
	//////////////particleFilter();
	//////////////// 调用 EMst 类中的其他函数
	//////////////g.CreateALGraph();
	/////////////////// ...


	/*step3:最小生成树*/
	GraphAdjList mindis;
	vector<int> unusedset;
	vector<int> cornerList;

	//构建一个等长数列
	int length = m_knnpoint->points.size();
	for (int i = 0; i < length; i++) cornerList.push_back(i);

	//邻接表赋予长度
	mindis.numVertsxes = length;
	//新建邻接表
	mindis.adjList = new AdjList[length];

	for (int i = 0; i < length; i++) {
		//初始化顶点分类
		mindis.adjList[i].data = cornerList[i];
		//初始化第一个边（之后建立邻接链表）
		mindis.adjList[i].firstedge = NULL;
		//初始化未使用列表
		unusedset.push_back(cornerList[i]);
	}
	//邻接表赋予
	mindis = g.Graph;

	//最小生成树的初始化，步骤同上
	GraphAdjList initial_tree;
	initial_tree.numVertsxes = length;
	initial_tree.adjList = new AdjList[length];
	for (int i = 0; i < initial_tree.numVertsxes; i++) {
		initial_tree.adjList[i].data = cornerList[i];
		initial_tree.adjList[i].firstedge = NULL;
		//设置度数为0
		initial_tree.adjList[i].degree = 0;
		//设置为未访问
		initial_tree.adjList[i].visit = false;
		//设置为非边角点
		initial_tree.adjList[i].is_cor = false;
	}

	/*Step3:Find the min span tree by Prim algothrim,but there is not only a  adjacency list*/
	//当未使用列表长度为空时结束循环,使用循环的原因是因为不止一个邻接表
	while (unusedset.size())
	{
		//初始化成本信息， vector.front()代表容器的第一个元素
		//这看上去效率可能有点问题，但是很方便
		int custar = order_element_in_vector(cornerList, unusedset.front());
		//获取初始化后的树地址
		GraphAdjList* tree1 = &initial_tree;
		//计算一遍最小生成树，结果存入了tree1（initial_tree）
		g.MiniTree_Prim(g.Graph, custar, tree1);
		// 计算度大于等于3的点有多少
		//初始化计数器
		int nodenum = 0;
		EdgeNode* v;
		//默认为不直
		bool is_notstraight = false;

		//这里的做法首先是记录有邻接的点，先存起来就行了
		vector<int> nodeNumList;

		//这个length要不要设置为final， 不过为什么是遍历长度为length呢，这样感觉每次都要重新扫描一遍
		for (int vec = 0; vec < length; vec++) {
			//这里把树的邻接表结点地址传递过来
			v = tree1->adjList[vec].firstedge;

			//这里从未使用表中擦除起点坐标
			unusedset.erase(std::remove(unusedset.begin(), unusedset.end(), cornerList[custar]), unusedset.end());
			//只要邻接表还在，就不会停止
			while (v) {
				//这里从未使用表中擦除邻接点的坐标
				unusedset.erase(std::remove(unusedset.begin(), unusedset.end(), tree1->adjList[v->adjvex].data), unusedset.end());
				//该点度+1
				tree1->adjList[vec].degree++;
				//下一个
				v = v->next;
			}
			if (tree1->adjList[vec].degree > 2) {
				nodenum++;
			}


		}
		//quickSort(nodeNumList, 0, nodeNumList.size() - 1);

		//KMST

			//这里有点奇怪，必须判断判断当前是否邻接存在
		for (int vec = 0; vec < length; vec++) {
			//这里把树的邻接表结点地址传递过来
			v = tree1->adjList[vec].firstedge;
			if (v)nodeNumList.push_back(vec);
		}
		/*
		首先选择起点，起点策略为选择第一个点（）；
		遍历该点找到最远的点（这里是图论的广度优先和深度优先）；
		以最远的点为起始点，
		找第二个最远的点，
		这两点之间的路径就是主干
		*/

		int A = searchByDistance(tree1, nodeNumList, 0);
		//如果A发生了错误就不进行搜索，-1表示未提前结束
		if (A == -1)
		{
			logger("本轮搜索失败");
		}
		else {
			//PathTree pt ;
			PathTree pt = searchByDistance(tree1, nodeNumList, A, 1);
			//策略 flag 0 只返回一个包含端点的结构体
			//flag 1 返回一个包含路径的结构体
			logger("本轮搜索完全" + s(pt.childChains.size()));
			outputSimplePolylineFile(g, m_knnpoint, pt, HeightLimit, Orderitreator);
		}
		logger("该轮结点总数为" + s(nodeNumList.size()));
		logger("该轮废弃结点总数为" + s(deprecatedList.size()));
		clear(nodeNumList);

	}
}



void txtToPCD(string txtname)
{
	int num_txt;
	cout << "Now is translating txt file to PCD and the point type is XYZ8D" << endl;
	//定义一种类型表示TXT中的点云格式xyz
	typedef struct TXT_Point_XYZ
	{
		double x;
		double y;
		double z;
		float intensity;
		float gvalue;
		float p0;
		float p1;
		float ptype;

	}TOPOINT_XYZ8D;

	//读取txt文件
	FILE *fp_txt;
	TXT_Point_XYZ txt_points;
	vector<TXT_Point_XYZ> my_vTxtPoints;
	fp_txt = fopen(txtname.c_str(), "r");

	if (fp_txt)
	{
		while (fscanf(fp_txt, "%lf %lf %lf %f %f %f %f %f", &txt_points.x, &txt_points.y, &txt_points.z, &txt_points.intensity,
			&txt_points.gvalue, &txt_points.p0, &txt_points.p1, &txt_points.ptype) != EOF)
			/*&txt_points.x, &txt_points.y, &txt_points.z,&txt_points.Tlambda2, &txt_points.Tlambda1,
				&txt_points.Tlambda0, &txt_points.gz, &txt_points.gy, &txt_points.gx, &txt_points.gvalue, &txt_points.pca_lamda2,
				&txt_points.pca_lamda1, &txt_points.pca_lamda0, &txt_points.intensity*/
		{//将点存入容器尾部
			my_vTxtPoints.push_back(txt_points);
		}
	}
	else cout << "读取txt文件失败" << endl;

	num_txt = my_vTxtPoints.size();

	//写入点云数据
	pcl::PointCloud<MyPointType> ::Ptr Bcloud(new pcl::PointCloud<MyPointType>);
	Bcloud->width = num_txt;
	Bcloud->height = 1;
	Bcloud->is_dense = false;
	Bcloud->points.resize(Bcloud->width * Bcloud->height);
	for (int i = 0; i < Bcloud->points.size(); ++i)
	{
		Bcloud->points[i].x = my_vTxtPoints[i].x;
		Bcloud->points[i].y = my_vTxtPoints[i].y;
		Bcloud->points[i].z = my_vTxtPoints[i].z;
		Bcloud->points[i].intensity = my_vTxtPoints[i].intensity;
		Bcloud->points[i].gvalue = my_vTxtPoints[i].gvalue;
		Bcloud->points[i].p0 = my_vTxtPoints[i].p0;
		Bcloud->points[i].p1 = my_vTxtPoints[i].p1;
		Bcloud->points[i].ptype = my_vTxtPoints[i].ptype;
	}
	string tname = txtname.substr(txtname.length() - 10, txtname.length() - 4);
	string suffix = ".pcd";
	string pcdname = tname + suffix;
	pcl::io::savePCDFileASCII("1.pcd", *Bcloud);
	cout << "从 txt_pcd.txt读取" << Bcloud->points.size() << "点写入txtTO.pcd" << endl;

}