#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include "ObjectiveRating.h"

using namespace std;

bool dataMatrixNormalized(vector<vector<double>>& indicesData, vector<double>& refeData, 
	vector<double>& EMS, vector<int> dataType) 
{
	if (dataType.size() != indicesData.size())
		return false;
	double lowBound = 0,upBound;
	for (int i = 0; i < indicesData.size(); i++) {
		auto upBoundpos = max_element(indicesData[i].begin(), indicesData[i].end());//返回相对于indices[i].begin()的位置
		upBound = (indicesData[i][upBoundpos - indicesData[i].begin()] > EMS[i] * 2 ? 
			indicesData[i][upBoundpos - indicesData[i].begin()] : EMS[i] * 2);//上边界为均方差2倍和最大值的较大者
		for (int j = 0; j < indicesData[i].size(); j++) {//正向指标和负向指标的标准化
			indicesData[i][j] = abs(abs((indicesData[i][j] - refeData[i])-lowBound) / (upBound - lowBound - dataType[i]));
		}
	}
	return true;
}

vector<double> EntropyWeightMethod(vector<vector<double>> indicesData) {
	vector <double> weight(indicesData.size(), 0);
	double sum = 0, sum_log = 0;
	for (int i = 0; i < indicesData.size(); i++) {
		sum = accumulate(indicesData[i].begin(), indicesData[i].end(), 0.0);
		for (int j = 0; j < indicesData[i].size(); j++) {
			if (indicesData[i][j]!=0)
				sum_log += -(indicesData[i][j] / sum) * log(indicesData[i][j] / sum);
		}
		weight[i] = -sum_log / log(indicesData[i].size());
		sum = 0;
		sum_log = 0;
	}
	sum = accumulate(weight.begin(), weight.end(), 0.0);
	for (vector<double>::size_type it = 0; it < weight.size(); it++) {
		weight[it] = sum + 1 - 2 * weight[it];
	}
	sum_log = accumulate(weight.begin(), weight.end(), 0.0);
	for (vector<double>::size_type it = 0; it < weight.size(); it++) {
		weight[it] = weight[it]/sum_log;
	}
	return weight;
}