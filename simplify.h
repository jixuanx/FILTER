#pragma once

#include "helper.h"
#include "EMst.h"

using namespace std;
//这里没有使用指针容器和链表，是因为考虑到如果给予地址，函数体运行后会丢失容器长度和映射等信息,不如使用递归方法可以避免这些问题
//这里特别需要注意，容器的push_back相当于是深拷贝了一个变量进去
typedef struct PathTree {
	vector<int> self;
	vector<PathTree> childChains;
}PathTree;

vector<int> deprecatedList;
//主干阈值
int PathMinLimit = 30;
//树大小阈值
int TreeMinLimit = 30;
//简化阈值
float HeightLimit = 0.005f;

//vector<int> DPSimplify(p_Mypoint &P, vector<int> &nodeNumList, vector<int> singlePath, float &tolerance);
//void outputSimplePolylineFile(p_Mypoint &P, PathTree &pt, float tolerance);
//
//template <typename PointT>
//float getHeight(PointT ps, PointT pe, PointT p0);
//
//ostream & operator<<(ostream &cout, PathTree &p);
//int searchByDistance(GraphAdjList* tree, vector<int> nodeNumList, int startPos);
//PathTree searchByDistance(GraphAdjList* tree, vector<int> nodeNumList, int startPos, int flag, p_Mypoint &P);
//p0表示求高的点

template <typename PointT>
float getHeight(PointT ps, PointT pe, PointT p0)
{
	float line_S0 = sqrt(pow(ps.x - p0.x, 2) + pow(ps.y - p0.y, 2) + pow(ps.z - p0.z, 2));
	float line_E0 = sqrt(pow(pe.x - p0.x, 2) + pow(pe.y - p0.y, 2) + pow(pe.z - p0.z, 2));
	float line_SE = sqrt(pow(ps.x - pe.x, 2) + pow(ps.y - pe.y, 2) + pow(ps.z - pe.z, 2));
	if (line_S0 + line_E0 > line_SE && line_S0 + line_SE > line_E0 && line_SE + line_E0 > line_S0) {
		{
			//计算周长
			float s = (line_S0 + line_E0 + line_SE) / 2;
			//海伦公式
			float area = sqrt(s*(s - line_S0)*(s - line_E0)*(s - line_SE));
			//计算p0对应高
			return area / line_SE * 2;
		}
	}
	else return min(line_S0, line_E0);
}

//输入：全局列表，单段列表，距离阈值
//辅助输入：树（储存坐标），结点编号表（映射使用）
//返回：简化路径
vector<int> DPSimplify(p_Mypoint &P, vector<int> &nodeNumList, vector<int> singlePath, float &tolerance)
{
	if (singlePath.size() > 2)
	{
		int ps = nodeNumList[singlePath[0]];
		int pe = nodeNumList[singlePath[singlePath.size() - 1]];

		int indexMax = 0;
		float hMax = 0;
		for (int i = 1; i < singlePath.size() - 1; i++)
		{
			int pi = nodeNumList[singlePath[i]];
			float dH = getHeight(P->points[ps], P->points[pe], P->points[pi]);
			if (dH > tolerance && dH > hMax)
			{
				hMax = dH;
				indexMax = i;
			}
		}
		//如果为0则中途点全部舍去
		if (hMax == 0)
		{
			singlePath.erase(singlePath.begin() + 1, singlePath.end() - 1);
			return singlePath;
		}
		else {
			vector<int> left;

			for (int i = 0; i < indexMax + 1; i++) {
				left.push_back(singlePath[i]);
			}

			vector<int> right;
			for (int i = indexMax; i < singlePath.size(); i++) {
				right.push_back(singlePath[i]);
			}
			vector<int> all;
			connect(all, DPSimplify(P, nodeNumList, left, tolerance));
			connect(all, DPSimplify(P, nodeNumList, right, tolerance));

			return all;
		}
	}
	else return singlePath;

}

template <typename PointT>
float getArea(PointT p1, PointT p2, PointT p3)
{
	float line_13 = sqrt(pow(p1.x - p3.x, 2) + pow(p1.y - p3.y, 2) + pow(p1.z - p3.z, 2));
	float line_23 = sqrt(pow(p2.x - p3.x, 2) + pow(p2.y - p3.y, 2) + pow(p2.z - p3.z, 2));
	float line_12 = sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));

	//计算周长
	float s = (line_13 + line_23 + line_12) / 2;
	//海伦公式
	float area = sqrt(s*(s - line_13)*(s - line_23)*(s - line_12));
	//计算p3对应高
	return area;

}

//输入：全局列表，单段列表，面积阈值
//辅助输入：树（储存坐标），结点编号表（映射使用）
//返回：简化路径
//输入：全局列表，单段列表，面积阈值
//辅助输入：树（储存坐标），结点编号表（映射使用）
//返回：简化路径
vector<int> VMSimplify(p_Mypoint &P, vector<int> &nodeNumList, vector<int> singlePath, float &tolerance)
{
	//这里不同的是因为使用了三角形，必须保证路径节点数大于3
	if (singlePath.size() > 3)
	{
		//singlePath是原始结点映射表，存储的是路径上的点在结点编号表中的位置，由于非引用传参，所以可以随意修改数据
		//对于三角形，按singlePath的顺序存储即可

		//初始化三角形面积顺序表，需要确保长度为结点映射表的长度-2
		vector<float> areaList;
		for (int i = 1; i < singlePath.size() - 1; i++) {
			int ps = nodeNumList[singlePath[i - 1]];
			int pm = nodeNumList[singlePath[i]];
			int pe = nodeNumList[singlePath[i + 1]];
			areaList.push_back(getArea(P->points[ps], P->points[pm], P->points[pe]));
		}

		//按顺序排列求出最小三角形所在位置和对应顶点
		//理论上来说这里由于只有局部发生了变化应该做动态更新，但是时间有限所以直接使用容器求最小值的函数.
		int count = 0;
		int baseSize = singlePath.size();
		while (true)
		{
			//求解最小三角形的下标
			int min_data = min_element(areaList.begin(), areaList.end()) - areaList.begin();
			//如果符合要求终止循环，或者出现异常循环时也需要终止
			if (areaList[min_data] > tolerance || singlePath.size() < 4)break;
			if (count > baseSize)break;
			//求解最小三角形对应顶点在映射表中的下标
			int pmin = min_data + 1;
			logger(singlePath[singlePath.size() - 1]);
			//logger("pos:",min_data,",size:",areaList.size(),",size2:", singlePath.size());
			//这里需要排除点位于第二个和倒数第二个的特殊情形，最后是位于中间的情形
			if ((pmin - 1) == 0)
			{
				int p1 = nodeNumList[singlePath[0]];
				//int pm = nodeNumList[singlePath[pmin]];
				int p3 = nodeNumList[singlePath[2]];
				int p4 = nodeNumList[singlePath[3]];
				float a1 = getArea(P->points[p1], P->points[p3], P->points[p4]);
				//更新面积列表
				areaList.erase(areaList.begin(), areaList.begin() + 2);
				areaList.insert(areaList.begin(), a1);
				//更新路径映射表
				singlePath.erase(std::begin(singlePath) + 1);
			}
			else if ((pmin + 1) == (singlePath.size() - 1))
			{
				int p1 = nodeNumList[singlePath[pmin - 2]];
				int p2 = nodeNumList[singlePath[pmin - 1]];
				//int pm = nodeNumList[singlePath[pmin]];
				int p4 = nodeNumList[singlePath[pmin + 1]];
				float a1 = getArea(P->points[p1], P->points[p2], P->points[p4]);
				//更新面积列表
				areaList.erase(areaList.end() - 2, areaList.end());
				areaList.push_back(a1);
				//更新路径映射表
				singlePath.erase(std::end(singlePath) - 2);
			}
			else {
				int p1 = nodeNumList[singlePath[pmin - 2]];
				int p2 = nodeNumList[singlePath[pmin - 1]];
				//int pm = nodeNumList[singlePath[pmin]];
				int p4 = nodeNumList[singlePath[pmin + 1]];
				int p5 = nodeNumList[singlePath[pmin + 2]];
				float a1 = getArea(P->points[p1], P->points[p2], P->points[p4]);
				float a2 = getArea(P->points[p2], P->points[p4], P->points[p5]);
				//更新面积列表
				areaList.erase(areaList.begin() + (min_data - 1), areaList.begin() + (min_data + 2));
				areaList.insert(areaList.begin() + (min_data - 1), a1);
				areaList.insert(areaList.begin() + (min_data), a2);
				//更新路径映射表
				singlePath.erase(singlePath.begin() + pmin);
			}
			if (singlePath.size() != areaList.size() + 2)logger(singlePath.size(), "//", areaList.size());
			count++;
		}


	}
	return singlePath;
}





template <typename PointT>
float getDistance(PointT p1, PointT p2)
{
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

//连续点阈值
float PointInterval = 0.01000f;
vector<MyPointType> createConsecutivePoints(const Mypoint &P, vector<int> &nodeNumList, vector<int> &pathList) {
	vector<MyPointType> pSet;
	int size = pathList.size();
	//压入起始点
	pSet.push_back(P[nodeNumList[pathList[0]]]);
	for (int i = 0; i < size-1; i++)
	{
		int pos1 = nodeNumList[pathList[i]];
		int pos2 = nodeNumList[pathList[i+1]];

		//生成一串连续点
		float d = getDistance(P[pos1], P[pos2]);
		//logger(d);
		//这里是中间点的数量
		int Num = round((d / PointInterval) - 1);
		float x0 = P[pos1].x >P[pos2].x ? P[pos1].x : P[pos2].x;
		float y0 = P[pos1].y > P[pos2].y ? P[pos1].y : P[pos2].y;
		float z0 = P[pos1].z > P[pos2].z ? P[pos1].z : P[pos2].z;
		if (Num > 0)
		{
			//重新分配间隔
			float dx = fabs((P[pos1].x - P[pos2].x) / ((float)Num));
			float dy = fabs((P[pos1].y - P[pos2].y) / ((float)Num));
			float dz = fabs((P[pos1].z - P[pos2].z) / ((float)Num));
			//logger(dx);
			for (int d=1; d <= Num; d++)
			{
				float x = x0 + d*dx;
				float y = y0 + d*dy;
				float z = z0 + d*dz;
				MyPointType dp;
				dp.x = x;
				dp.y = y;
				dp.z = z;
				pSet.push_back(dp);
			}
		}

		//将第二段终点压入
		pSet.push_back(P[pos2]);
	}

	return pSet;

}
void adjustProbability(EMst &g, vector<int> &nodeNumList, vector<int> &pathList)
{
	
	//logger(nodeNumList);
	//计算最小距离
	int size = nodeNumList.size();

	vector<MyPointType> pSet = createConsecutivePoints(g.knnPoint, nodeNumList, pathList);
	int childSize = pSet.size();
	for (int i = 0; i < size; i++)
	{
		int pos = nodeNumList[i];
		float baseD = Inf;
		
		for (int j = 0; j < childSize; j++)
		{
			float d = getDistance(g.knnPoint[pos], pSet[j]);
			if (d < baseD)baseD = d;
		}
		float p0 = 0.5 * (g.knnPoint[pos].p0 + exp(-(pow(baseD, 2))));
		//cout << "更新p0：" << p0 << endl;
		//if (p0 > threshold_p) 
		//{
		//	nodeNumList.front()
		//	g.knnPoint.erase(i); i--; }
		g.knnPoint[pos].p0 = p0;
		g.knnPoint[pos].p1 = 1 - p0;
		g.knnPoint[pos].ptype = p0 > 0.5 ? 0 : 1;
	}
}

//这个函数是为了查询元素在容器数组里的次数
int order(vector<int> v, int element) {
	for (int order = 0; order < v.size(); order++) {
		if (element == v[order])
			return order;
	}
	//为意外情况做个补充
	return -1;
}

void outputSimplePolylineFile(EMst &g, p_Mypoint &P, PathTree &pt, float tolerance,int Orderitreator)
{
	vector<int> all = pt.self;
	if (pt.childChains.size() > 0)
	{
		//这也是一个映射，储存的是path里面的的真实坐标序列编号
		vector<int> rootNodeList;

		for (int i = 0; i < pt.childChains.size(); i++)
		{

			vector<int> c = pt.childChains[i].self;
			int pos = order(all, c[c.size() - 1]);
			if (pos != -1) {
				rootNodeList.push_back(pos);
			}
			else {
				logger("异常" + s(all[0]));
			}
			outputSimplePolylineFile(g, P, pt.childChains[i], tolerance, Orderitreator);
		}
		//整理顺序，方便排列
		quickSort(rootNodeList, 0, rootNodeList.size() - 1);
		vector<int> singlePath;
		//outputPolylineFile(P, pt.self, createSequentialList(0, pt.self.size()), s(pt.self[0]));
		for (int i = 0; i < all.size(); i++)
		{
			singlePath.push_back(all[i]);
			if ((rootNodeList.size() > 0 && i == rootNodeList[0]) ||
				i == all.size() - 1)
			{
				PathTree pt0 = *new PathTree();
				pt0.self = singlePath;
				//pt0.root = vectorToList(singlePath);
				outputSimplePolylineFile(g, P, pt0, tolerance, Orderitreator);

				if (i < all.size() - 1) {
					rootNodeList.erase(rootNodeList.begin());
					singlePath.clear();
					singlePath.push_back(all[i]);
				}

			}
		}
	}
	else
	{
		vector<int> pathMap = createSequentialList(0, all.size());

		//debugList(pt.self);
		if (all.size() > 0)
		{
			//可以再传一个未简化版本过去简化计算
			//outputPolylineFile(P, pt.self, pathMap, s(pt.self[0]));
			//vector<int> resultPath = DPSimplify(P, all, pathMap, tolerance);
			vector<int> resultPath = VMSimplify(P, all, pathMap, tolerance);
			outputPolylineFile(P, all, resultPath, s(Orderitreator));
			adjustProbability(g, all, resultPath);
		}
	}
}


ostream & operator<<(ostream &cout, PathTree &p)
{
	cout << "[";
	int i = 0;
	for (; i < p.self.size() - 1; i++)cout << p.self[i] << ",";
	if (i == p.self.size() - 1)cout << p.self[i] << "]";
	if (p.childChains.size() > 0)
	{
		cout << "[";
		for (int j = 0; j < p.childChains.size(); j++) {
			PathTree pt = p.childChains[j];
			cout << pt << ",";
		}
		cout << "]";
	}
	return cout;
}

/*
这里选择广度优先的原因是因为无论如何，我们都要遍历所有的分支以避免漏
广度优先算法：

*/
//startPos这里应该指的是在node列表里面的位置
//nodeNumList存储对应位置的点在全图中的邻接表位置
//返回值为也是在nodeNumList的位置
void dfs(GraphAdjList* tree, vector<int> *nodeNumList, vector<int> *visitedList, int root, int dep, vector<int> *depthList) {
	int childDep = dep + 1;
	//cout << root << endl;
	EdgeNode *first = tree->adjList[nodeNumList->at(root)].firstedge;
	while (first)
	{
		if (order_element_in_vector(*(vector<int>*) visitedList, first->adjvex) == -1)
		{
			int childPositon = order_element_in_vector(*(vector<int>*) nodeNumList, first->adjvex);
			if (childPositon != -1)
			{
				visitedList->push_back(first->adjvex);
				dfs(tree, nodeNumList, visitedList, childPositon, childDep, depthList);
				//为什么这里可以写成下标0，即depthList[0][childPositon]
				(*depthList)[childPositon] = childDep;
				//cout << first->adjvex << "子节点" <<","<< childDep << depthList->at(childPositon) << endl;
			}
		}
		first = first->next;
	}
}
//这里找路径其实先记录深度，然后向上寻找深度低的地方即可，广度和深度优先都能表示
//这里深度回溯不确定要不要优化起始点选择（若出错改为单邻居的起始点）
//结果方向为最深度至起点
vector<int> backtracking(GraphAdjList* tree, vector<int> &nodeNumList, vector<int> &depthList, int maxDepthPos, int root) {

	//cout << nodeNumList[maxDepthPos] << endl;
	//深度
	int dep = depthList[maxDepthPos];
	//这里表示最深的结点的父结点
	//这里有个问题是没有把未考虑起点的其他邻接
	EdgeNode *endFather = tree->adjList[nodeNumList[maxDepthPos]].firstedge;
	while (endFather) {
		int nPos = order_element_in_vector(nodeNumList, endFather->adjvex);
		if (depthList[nPos] < depthList[maxDepthPos])
		{
			dep = depthList[nPos];

			break;
		}
		endFather = endFather->next;
	}

	//路径列表，由后往前
	vector<int> pathList;
	//将第一个点压入路径列表
	pathList.push_back(maxDepthPos);
	//这里的运行条件为大于0
	while (dep > 0)
	{
		//首先将该父节点压入路径列表
		pathList.push_back(order_element_in_vector(nodeNumList, endFather->adjvex));
		//这里是遍历父节点的邻居，假如父节点有一个深度比父节点小的点，则其为父节点的父节点
		vector<EdgeNode *> neighborList;
		//选择父节点的邻居
		EdgeNode *friendNode = tree->adjList[endFather->adjvex].firstedge;
		while (friendNode)
		{
			//记录下来邻居
			neighborList.push_back(friendNode);
			//前往下一个邻居
			friendNode = friendNode->next;
		}
		//寻找深度小的邻居
		for (int i = 0; i < neighborList.size(); i++)
		{
			int neighborPos = order_element_in_vector(nodeNumList, neighborList[i]->adjvex);
			if (depthList[neighborPos] < dep)
			{
				//order_element_in_vector(nodeNumList, friendNode->adjvex)
				dep = depthList[neighborPos];
				endFather = neighborList[i];
				break;
			}
		}
		//清理邻居表空间
		neighborList.clear();
	}
	//由于起点的深度问题，所以这里需要补充起点
	pathList.push_back(root);
	return pathList;
}




//这里flag表示是否输出路径
PathTree searchByDistance(GraphAdjList* tree, vector<int> nodeNumList, int startPos, int flag) {
	logger(s(flag));
	int totalAmount = nodeNumList.size();

	vector<int> visitedList;

	//depthList和nodeNumList是一个映射关系，同位元素是对某对象的不同属性的表达，这里后面也可以考虑用结构体处理，
	//因此禁止对容器中点的数量进行简单操作，如果要操作，需要进行联动操作
	vector<int> depthList;
	for (int vec = 0; vec < nodeNumList.size(); vec++)depthList.push_back(-1);

	visitedList.push_back(nodeNumList[startPos]);
	//tagList[nodeNumList[startPos]] = 1;
	depthList[startPos] = 0;
	dfs(tree, &nodeNumList, &visitedList, startPos, 0, &depthList);


	//假如只需要获取最远点，那么这里先确定
	int maxDepthPos = 0;
	for (int i = 0; i < depthList.size(); i++)
	{
		if (depthList[i] > depthList[maxDepthPos])maxDepthPos = i;
	}
	visitedList.swap(vector<int>());

	//if (flag == 0)return depthList.size() > 0 ? maxDepthPos : -1;
	if (flag == 0) {
		PathTree tmp0;
		vector<int> tmp; tmp.push_back(maxDepthPos);
		tmp0.self = tmp;
		return tmp0;
	}

	//注意路径表需要配合结点映射表使用
	vector<int> pathList = backtracking(tree, nodeNumList, depthList, maxDepthPos, startPos);
	depthList.swap(vector<int>());
	//如果这样的话做个改进，根据flag决定返回的链表
	//if (flag == 0)return *new vector<int>->push_back(maxDepthPos);

	//不管怎么样，都需要禁忌路径结果，区别在于如果路径过短则禁忌全部，注意路径表需要配合结点映射表使用
	PathTree childPath;
	vector<int> realPath = pathList;
	if (pathList.size() < PathMinLimit)
	{
		connect(deprecatedList, nodeNumList);
		//connectByMapping(deprecatedList, pathList, nodeNumList);

		return *new PathTree();
	}
	else
	{
		inverseMapping(realPath, nodeNumList);

		childPath.self = realPath;

		//debugList(realPath);
		connectByMapping(deprecatedList, pathList, nodeNumList);
	}


	//输出路径到文件中
	string filename = to_string(nodeNumList[startPos]);
	//outputPolylineFile(P, nodeNumList, pathList, filename);
	//logger(filename);

	/*存入pathList的点表示已经被标记为主干的点，但是这个列表并不能直接被我们使用，因为它的目的是作为前后点序输出存在
	我们需要从nodeNumList剩下的点中进行排查，这些点都保有深度，且深度与主干息息相关
	这些点存入新的结点编号表，
	主干上的结点就是分支的起点，也是分支深度最低的地方
	*/

	//只有当flag为2时才执行分支搜索
	/*
	如果要处理分支的分支这种情况，那么必须进行筛选，因此现在开始讨论，这里给出两个先行限定约束
	第一个，子树的根节点是其主树的主干结点，是子树主干的起点；
	第二个，有两种情况的节点不需要被考虑：
		第一种是已经查明其所在子树点太少，那么该子树上所有点置入废弃集合（若最终输出点数过少同样全部放入该集合）【废弃集合可以考虑从全结点映射编号表上进行移除，但是还需要同时考虑在映射表中查询点位置时返回为-1情况】
		第二种是已经验证处于主干的点（注意这种情况需要特别保留结点表）。

	子树查询的过程是一个递归过程，大致如下
	1级主干查询结点，随后编制子树映射编号表，若表太小则禁忌，否则探索子树最深处节点，回溯查询深度是否达标，达标则输出路径，否则全部禁忌；
		n级主干现在有当前主干的节点顺序表，同第一步处理，若无足够大的子树则结束在n级主干的查询
	强制终止条件：没有点在未归类集合


	*/
	if (flag > 0)
		for (int i = 0; i < pathList.size(); i++)
		{
			int adjvexMainBranch = nodeNumList[pathList[i]];
			//度大于2才可能是强壮分支的起点
			if (tree->adjList[adjvexMainBranch].degree > 2)
			{

				//查询起点邻居
				EdgeNode *friendNode = tree->adjList[adjvexMainBranch].firstedge;
				while (friendNode)
				{

					//对于子节点邻居而言，由于前进深度可能不变（如果在同一层子子结点存在多个），这里需要再搜索一次最深的子节点
					//一个暂时性的解决方法，我们需要提供子树的结点编号表，由于移除了主干，需要重新构树，这里可以直接用广度优先处理
					int friendAddressGlobal = friendNode->adjvex;
					//只有当邻居不在主干时才能作为子节点
					if (order_element_in_vector(pathList, order_element_in_vector(nodeNumList, friendAddressGlobal)) == -1
						&& order_element_in_vector(deprecatedList, friendAddressGlobal) == -1)
					{

						queue<int> nodeBranchQueque;
						//这里在检索生成子树的结点编号表，因此我们需要使用全局下标，而不是结点标号表的映射下标
						nodeBranchQueque.push(friendAddressGlobal);
						vector<int> childNodeNumList;
						childNodeNumList.push_back(friendAddressGlobal);
						while (nodeBranchQueque.size()) {

							//探索子节点
							queue<int> childLayerQueque;
							while (nodeBranchQueque.size())
							{
								EdgeNode *first = tree->adjList[nodeBranchQueque.front()].firstedge;
								while (first)
								{
									if (first->adjvex != adjvexMainBranch
										&& order_element_in_vector(childNodeNumList, first->adjvex) == -1
										)
									{
										childLayerQueque.push(first->adjvex);
										childNodeNumList.push_back(first->adjvex);
									}
									first = first->next;
								}
								nodeBranchQueque.pop();
							}
							nodeBranchQueque = childLayerQueque;
							childLayerQueque.empty();
						}
						//这里最好设置为全局阈值
						if (childNodeNumList.size() > TreeMinLimit)
						{	//理论上应该是0，是为了避免出错，查一下位置（由于是顺序找，会很快）
							childNodeNumList.insert(childNodeNumList.begin(), adjvexMainBranch);

							PathTree reslut = searchByDistance(tree, childNodeNumList, order_element_in_vector(childNodeNumList, adjvexMainBranch), 2);
							if (reslut.self.size() > 0)
								childPath.childChains.push_back(reslut);

							childNodeNumList.swap(vector<int>());
							//cout << childNodeNumList[B] << endl;
						}
						else
						{
							connect(deprecatedList, childNodeNumList);
						}
					}
					//前往下一个邻居
					friendNode = friendNode->next;
				}
			}
		}
	if (flag > 0)
	{

		return childPath;
	}
	return  *new PathTree();

}

int searchByDistance(GraphAdjList* tree, vector<int> nodeNumList, int startPos) {
	vector<int> result = searchByDistance(tree, nodeNumList, startPos, 0).self;
	return result.size() > 0 ? result[0] : -1;
}