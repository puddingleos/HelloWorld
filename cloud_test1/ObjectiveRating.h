#pragma once
#ifndef _OBJECTIVERATING_H_
#define _OBJECTIVERATING_H_

#include <vector>

using namespace std;

bool dataMatrixNormalized(vector<vector<double>>& indicesData, vector<double>& refeData, vector<double>& EMS, vector<int> dataType);
vector<double> EntropyWeightMethod(vector<vector<double>> indicesData);
vector<double> GreyRelationalAnalysis(vector<vector<double>> indicesData);
vector<double> PrincipleComponentAnalysis(vector<vector<double>> indicesData);
#endif
