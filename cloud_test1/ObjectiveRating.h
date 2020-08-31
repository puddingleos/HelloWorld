#pragma once
#ifndef _OBJECTIVERATING_H_
#define _OBJECTIVERATING_H_

#include <vector>

using namespace Eigen;
using namespace std;

bool dataMatrixNormalized(vector<vector<double>>& indicesData, vector<double>& refeData, vector<double>& EMS, vector<int> dataType);
// Objective Weighting Method
vector<double> EntropyWeightMethod(vector<vector<double>> indicesData);
vector<double> GreyRelationalAnalysis(vector<vector<double>> indicesData);
vector<double> PrincipleComponentAnalysis(vector<vector<double>> indicesData);
vector<double> EvlWithTOPSIS(vector<double> weight, vector<vector<double>> normIndicesData);
//Plot
void PlotWithDataLength(vector<vector<double>> weight, vector<vector<int>> dataLength, 
	vector<vector<vector<string>>> IndexName, const string title = {});
void PlotWithEMS(vector<vector<double>> wEWM, vector<vector<double>> EMS, vector<vector<int>>  dataLength,
	vector<vector<vector<string>>> IndexName, const string title = {});

vector<vector<vector<double>>> dataConvert(vector<vector<double>> indicesData, string subRateM);
vector<double> AHP(vector<vector<vector<double>>> indicesData, vector<vector<string>>& indexName);
vector<double> consistencyCheck(vector<double> indicesMatrix, vector<vector<string>>& indexName);
vector<double> FDM(vector<vector<vector<double>>> indicesData, vector<vector<string>>& indexName);

bool PyCvxpyInputData(
	//主观评价权重和数据
	string submethods,
	vector<double> subweight,
	vector<vector<double>> subindicesData,
	//客观评价权重和数据
	string objmethods,
	vector<double> objweight,
	vector<vector<double>> objindicesData,
	//返回结果
	vector<vector<double>>& A,
	vector<double>& b);

void polyfit_test(int n, double x[], double y[], int poly_n, double p[]);
void gauss_solve(int n, double A[], double x[], double b[]);

#endif
