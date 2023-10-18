#pragma once

#include <stdio.h>
//#include <stdarg.h>
//#include <direct.h> 
#include "EMst.h"
#include<tuple>

//template<typename T>
//void debug_notime(T t);
//template<typename T, typename... U>
//void debug_notime(T t, U... ts);
//template<typename T, typename... U>
//void logger(T t, U... ts);
//
//string s(int s);
//string  ss(int num);
//void quickSort(vector<int> &arr, int _left, int _right);
//void outputLineFile(p_Mypoint &P, vector<int> &nodeNumList, vector<int> &pathList, string name);
//void outputPolylineFile(p_Mypoint &P, vector<int> &nodeNumList, vector<int> &pathList, string name);
//void cleanFile(string strPath);
//void clear(vector<int> &v);
//void connect(vector<int> &v1, vector<int> &v2);
//void connectByMapping(vector<int> &v1, vector<int> &v2, vector<int> &map);
//void inverseMapping(vector<int> &v1, vector<int> &map);
//vector<int> createSequentialList(int first, int length);
//ostream & operator<<(ostream &cout, vector<int> &pathList);

template<typename T>
void clear(vector<T> &v)
{
	v.swap(vector<T>());
}

string s(int s)
{
	return to_string(s);
}

string ss(int num)
{
	return num > 9 ? s(num) : s(0) + s(num);
}


ostream & operator<<(ostream &cout, vector<int> &list) {
	//测试结果
	string x = "";
	for (int i = 0; i < list.size(); i++)
	{
		x = x + to_string(list[i]) + ",";
	}
	cout << x;
	return cout;
}


template<typename T>
void debug_notime(T t)
{
	std::cout << t << endl;
}

template<typename T, typename... U>
void debug_notime(T t, U... ts) {
	std::cout << t ;
	debug_notime(ts...);
}

template<typename T, typename... U>
void logger(T t, U... ts) 
{
	time_t time_seconds = time(0);
	struct tm t_tm;
	localtime_s(&t_tm, &time_seconds);
	string   str = "[" + ss(t_tm.tm_hour) + ":" + ss(t_tm.tm_min) + ":" + ss(t_tm.tm_sec) + "] ";
	std::cout << str;
	debug_notime(t, ts...);
}

void quickSort(vector<int> &arr, int _left, int _right) {
	int left = _left;
	int right = _right;
	int temp = 0;
	if (left < right) {   //待排序的元素至少有两个的情况
		temp = arr[left];  //待排序的第一个元素作为基准元素
		while (left != right) {   //从左右两边交替扫描，直到left = right

			while (right > left&& arr[right] >= temp)
				right--;        //从右往左扫描，找到第一个比基准元素小的元素
			arr[left] = arr[right];  //找到这种元素arr[right]后与arr[left]交换

			while (left < right && arr[left] <= temp)
				left++;         //从左往右扫描，找到第一个比基准元素大的元素
			arr[right] = arr[left];  //找到这种元素arr[left]后，与arr[right]交换

		}
		arr[right] = temp;    //基准元素归位
		quickSort(arr, _left, left - 1);  //对基准元素左边的元素进行递归排序
		quickSort(arr, right + 1, _right);  //对基准元素右边的进行递归排序
	}
}

void connectByMapping(vector<int> &v1, vector<int> &v2, vector<int> &map)
{
	for each (int num in v2)
	{
		v1.push_back(map[num]);
	}
}

void inverseMapping(vector<int> &v1, vector<int> &map) {
	for (int i = 0; i < v1.size(); i++)
	{
		v1[i] = map[v1[i]];
	}

}

void connect(vector<int> &v1, vector<int> &v2)
{
	v1.insert(v1.end(), v2.begin(), v2.end());
}



void outputLineFile(p_Mypoint &P, vector<int> &nodeNumList, vector<int> &pathList, string name)
{
	string filename = path;
	ofstream PLfile1(filename, ios::app);
	for (int i = 0; i < pathList.size() - 1; i++) {
		int start = nodeNumList[pathList[i]];
		int end = nodeNumList[pathList[i + 1]];
		PLfile1 << "GOCAD PLine 1\n";
		PLfile1 << "HEADER{\n";
		PLfile1 << "name:" << "line-" << name << "\n";
		PLfile1 << "}\n";
		PLfile1 << "ILINE\n";
		PLfile1 << "VRTX " << 0 << " " << P->points[start].x << " " << P->points[start].y << " " << P->points[start].z << "\n";
		PLfile1 << "VRTX " << 1 << " " << P->points[end].x << " " << P->points[end].y << " " << P->points[end].z << "\n";
		PLfile1 << "SEG " << 0 << " " << 1 << "\n";
		PLfile1 << "END";
	}
	PLfile1.close();
};

void outputPolylineFile(p_Mypoint &P, vector<int> &nodeNumList, vector<int> &pathList, string name)
{
	string filename = path + name +"3c04-30-30-vw0.005.pl";
	ofstream PLfile1(filename, ios::app);
	PLfile1 << "GOCAD PLine 0.02sim\n";
	PLfile1 << "HEADER{\n";
	PLfile1 << "name:" << "als-" << name << "\n";
	PLfile1 << "}\n";
	PLfile1 << "ILINE\n";
	for (int i = 0; i < pathList.size(); i++) {
		int start = nodeNumList[pathList[i]];
		PLfile1 << "VRTX " << i << " " << P->points[start].x << " " << P->points[start].y << " " << P->points[start].z << "\n";
	}
	for (int i = 0; i < pathList.size() - 1; i++) {

		PLfile1 << "SEG " << i << " " << i + 1 << "\n";

	}
	//由于最终输出路径都会经过这里，这里可以计算出来中间点，在树里分析（特别注意的是，这里有结点编号表），
	//就可以计算每个点对应的编号了。
	//在这里更新P0是完全可行的
	
	PLfile1 << "END";
	PLfile1.close();
};


void cleanFile(string strPath) 
{
	ofstream PLfile1(strPath);
	PLfile1.close();
}

vector<int> createSequentialList(int first, int length)
{
	vector <int> v;
	for (int i = 0; i < length; i++)
	{
		v.push_back(first + i);
	}
	return v;
}

//这个函数是为了查询元素在容器数组里的次数
int order_element_in_vector(const vector<int> &v, int element) {
	for (int order = 0; order < v.size(); order++) {
		if (element == v[order])
			return order;
	}
	//为意外情况做个补充
	return -1;
}
