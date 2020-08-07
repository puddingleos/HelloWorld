#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include "ObjectiveRating.h"

//#define DEBUG

using namespace std;
using namespace Eigen;



//��׼������
bool dataMatrixNormalized(vector<vector<double>>& indicesData, vector<double>& refeData, 
	vector<double>& EMS, vector<int> dataType) 
{
	if (dataType.size() != indicesData.size())
		return false;
	double lowBound = 0,upBound;
	for (int i = 0; i < indicesData.size(); i++) {
		auto upBoundpos = max_element(indicesData[i].begin(), indicesData[i].end());//���������indices[i].begin()��λ��
		upBound = (indicesData[i][upBoundpos - indicesData[i].begin()] > EMS[i] * 2 ? 
			indicesData[i][upBoundpos - indicesData[i].begin()] : EMS[i] * 2);//�ϱ߽�Ϊ������2�������ֵ�Ľϴ���
		for (int j = 0; j < indicesData[i].size(); j++) {//����ָ��͸���ָ��ı�׼��
			indicesData[i][j] = abs((abs(indicesData[i][j] - refeData[i])-lowBound) / (upBound - lowBound) - dataType[i]);
		}
	}
	return true;
}

//��Ȩ����ȨEWM
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

//�Ҷȹ�����GRA
vector<double> GreyRelationalAnalysis(vector<vector<double>> indicesData) {
	vector <double> weight(indicesData.size(), 0);
	vector<double> maxValue, minValue;
	for (vector<double>::size_type it = 0; it < indicesData.size(); it++) {
		for (vector<double>::size_type jt = 0; jt < indicesData[it].size(); jt++) {
			indicesData[it][jt] -= indicesData[it][indicesData[it].size()-1];// �ο�ֵ�������۾��󣨲ο�ֵһ��ѡ��һ�л����һ�У�
			if (indicesData[it][jt] < 0)
				indicesData[it][jt] = -indicesData[it][jt];
		}
		maxValue.push_back(indicesData[it][max_element(indicesData[it].begin(), indicesData[it].end()) - indicesData[it].begin()]);
		minValue.push_back(indicesData[it][min_element(indicesData[it].begin(), indicesData[it].end()) - indicesData[it].begin()]);
	}
	double maxMaxValue = maxValue[max_element(maxValue.begin(), maxValue.end()) - maxValue.begin()];
	double minMinValue = minValue[min_element(maxValue.begin(), maxValue.end()) - maxValue.begin()];
	for (vector<double>::size_type it = 0; it < weight.size(); it++) {
		for (vector<double>::size_type jt = 0; jt < indicesData[it].size(); jt++)
			weight[it] += (minMinValue + maxMaxValue * 0.5) / (indicesData[it][jt] + maxMaxValue * 0.5);
		weight[it] /= indicesData[it].size();
	}
	double sum = accumulate(weight.begin(), weight.end(), 0.0);
	for (vector<double>::size_type it = 0; it < weight.size(); it++)
		weight[it] /= sum;
	return weight;
}

// ���ɷַ�����
vector<double> PrincipleComponentAnalysis(vector<vector<double>> indicesData) {
	vector<double> weight(indicesData.size(), 0);
	vector<int> dataLength;
	vector<vector<double>> indicesData_t;
	// ����С����������ȡ����
	for (vector<double>::size_type it = 0; it < indicesData.size(); it++)
		dataLength.push_back(indicesData[it].size());
	int minLength = dataLength[min_element(dataLength.begin(), dataLength.end()) - dataLength.begin()];

	//����Eigen�⺯�����о�������
	MatrixXd indicesMatrix((int)dataLength.size(), minLength);
	for (vector<double>::size_type it = 0; it < indicesData.size(); it++) {
		for (vector<double>::size_type jt = 0; jt < minLength; jt++) {
			indicesMatrix(it, jt) = indicesData[it][jt];
		}
	}
	//1. denoisng
	//MatrixXd meanData = indicesMatrix.rowwise().mean();
	MatrixXd deMeanMatrix = indicesMatrix - indicesMatrix.rowwise().mean() * MatrixXd::Ones(1, minLength);
	indicesMatrix = deMeanMatrix.array() / (deMeanMatrix.rowwise().norm()* MatrixXd::Ones(1, minLength)).array() /(double)minLength;
	//2. Э����
	MatrixXd CovMatrix;
	CovMatrix = indicesMatrix * indicesMatrix.transpose();
#ifdef DEBUG
	cout << CovMatrix << endl;
#endif
	//3. �����ֽ�
	EigenSolver<MatrixXd> es;
	es.compute(CovMatrix, true);
	
	MatrixXi sortList = MatrixXi::Zero(es.eigenvalues().size(), 1);
	MatrixXi visited = MatrixXi::Zero(es.eigenvalues().size(),1);
	MatrixXd eigvalueSortDescend = MatrixXd::Zero(es.eigenvalues().size(),1);
	for (int i = 0; i < es.eigenvalues().size(); i++) {
		int it = 0;
		while (visited(it)) {
			if (it >= es.eigenvalues().size())
				it = 0;
			else
				it++;
		}
		eigvalueSortDescend(i) = abs(es.eigenvalues()(it));
		for (int j = 0; j < es.eigenvalues().size(); j++) {
			if (eigvalueSortDescend(i) <= abs(es.eigenvalues()(j)) && visited(j)!=1){
				eigvalueSortDescend(i) = abs(es.eigenvalues()(j));
				sortList(i) = j;
			}
		}
		visited(sortList(i)) = 1;
	}
#ifdef DEBUG
	cout << es.eigenvalues().cwiseAbs() << endl;
	cout << sortList << endl;
	cout << eigvalueSortDescend << endl;
#endif
	//4. ������ɸѡ
	int i = 0;
	double sum_t = 0;
	while (sum_t / eigvalueSortDescend.sum() < 0.85 && i < eigvalueSortDescend.size())
		sum_t += eigvalueSortDescend(i++);
	//5. �����¿ռ�
	MatrixXcd eigVectorNew = MatrixXcd::Zero(es.eigenvalues().size(), i);
	MatrixXd indicesDataSelected = MatrixXd::Zero(es.eigenvalues().size(), i);
	for (int j = 0; j < i; j++) {
		eigVectorNew.col(j) << es.eigenvectors().col(sortList(j));
		indicesDataSelected.col(j) << CovMatrix.col(sortList(j));
	}
	//cout << "eigVectorNew: " << eigVectorNew.transpose() << endl;
	//cout << "indicesDataSelected: " << indicesDataSelected << endl;
	MatrixXd T = (eigVectorNew.transpose() * indicesDataSelected).cwiseAbs2().real();

	////6. ����ָ��Ȩ��
	MatrixXd Ttmp = (T - T.rowwise().mean() * MatrixXd::Ones(1, T.rows())).rowwise().norm();
	Ttmp = Ttmp.array() / (Ttmp.sum() * MatrixXd::Ones(Ttmp.rows(), 1)).array();

	//cout << Ttmp << endl;

	for (int i = 0; i < Ttmp.rows(); i++) {
		weight[sortList(i)] = Ttmp(i);
	}

	return weight;
}