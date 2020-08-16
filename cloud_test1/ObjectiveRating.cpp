#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

#include "matplotlibcpp.h"
#include "ObjectiveRating.h"
#include <regex>

//#define DEBUG

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;


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


vector<double> EvlWithTOPSIS(vector<double> weight, vector<vector<double>> normIndicesData) {
	vector<int> dataLength;
	vector<double> EvalueResult;
	MatrixXd weight_t(weight.size(), 1);

	// ����С����������ȡ����
	for (vector<double>::size_type it = 0; it < normIndicesData.size(); it++)
		dataLength.push_back(normIndicesData[it].size());
	int minLength = dataLength[min_element(dataLength.begin(), dataLength.end()) - dataLength.begin()];
	//����Eigen�⺯�����о�������
	MatrixXd indicesMatrix((int)dataLength.size(), minLength);
	for (vector<double>::size_type it = 0; it < normIndicesData.size(); it++) {
		weight_t(it, 0) = weight[it];
		for (vector<double>::size_type jt = 0; jt < minLength; jt++) {
			indicesMatrix(it, jt) = normIndicesData[it][jt];
		}
	}
	//1. ��׼����Ȩ
	indicesMatrix = indicesMatrix.array() * (weight_t * MatrixXd::Ones(1, minLength)).array();
	//2. ������������
	MatrixXd PostiveIdeaSolution = indicesMatrix.rowwise().maxCoeff();
	MatrixXd NegativeIdeaSolution = indicesMatrix.rowwise().minCoeff();
	MatrixXd Dp = (indicesMatrix - PostiveIdeaSolution * MatrixXd::Ones(1, minLength)).colwise().norm();//�������֮��
	MatrixXd Dn = (indicesMatrix - NegativeIdeaSolution * MatrixXd::Ones(1, minLength)).colwise().norm();
	//3. ��������������(�������=1-�������)
	MatrixXd Cp = Dp.array() / (Dp + Dn).array();

	for (int i = 0; i < Cp.cols(); i++)
		EvalueResult.push_back(Cp(i));

	return EvalueResult;
}






void PlotWithDataLength(vector<vector<double>> wEWM, vector<vector<int>> dataLength, 
	vector<vector<vector<string>>> IndexName, const string title) {
	string IndexName_t;
	vector<double> wEWM_t;
	vector<double> x_t;
	VectorXd weight_t((int)wEWM.size());//����Eigen���������
	double meanValue(0.0), normValue(0.0), minValue(0.0), maxValue(0.0);

	int plotCount10 = 0;
	string plotType = "*--";
	for (int j = 0; j < wEWM[0].size(); j++) {//ָ���������0<j<=12
		for (int i = 0; i < wEWM.size(); i++) {//�ļ�������������������� 0<i<10
			wEWM_t.push_back(wEWM[i][j]);
			weight_t(i) = wEWM[i][j];
			x_t.push_back((double)dataLength[i][0]);//ѡȡ��һ��ָ���������Ϊ�ο�
		}
		for (int k = 0; k < IndexName[0][j].size(); k++)
			IndexName_t.append(IndexName[0][j][k]);

		//���ݷ���
		meanValue = weight_t.mean();//ƽ��ֵ
		normValue = (weight_t - meanValue * VectorXd::Ones(weight_t.size())).norm();//��׼��
		minValue = weight_t.minCoeff();//��Сֵ
		maxValue = weight_t.maxCoeff();//���ֵ

		IndexName_t.append(":"+to_string(meanValue)+","+to_string(normValue));

		if (plotCount10++ >= 10)
			plotType[0] = 'o';
		plt::named_plot(IndexName_t.c_str(), x_t, wEWM_t, plotType);

		IndexName_t.clear();
		wEWM_t.clear();
		x_t.clear();
	}
	plt::title(title);
	plt::xlabel("Amounts of data");
	plt::ylabel("Weight");
	plt::legend();
	plt::show();
}


void PlotWithEMS(vector<vector<double>> wEWM, vector<vector<double>> EMS, vector<vector<int>>  dataLength,
	vector<vector<vector<string>>> IndexName, const string title) {
	string dataLabel, xticks_t;
	vector<double> wEWM_t;
	vector<double> x_label, x_t, x_sort;
	vector<int> visited, sortList;
	vector<string> x_ticks;
	VectorXd weight_t((int)wEWM.size());//����Eigen���������
	double meanValue(0.0), normValue(0.0), minValue(0.0), maxValue(0.0);

	int plotCount10 = 0;
	string plotType = "*--";

		
	for (int i = 0; i < wEWM.size(); i++) {//�ļ�������������������� 0<i<10
		x_label.clear();
		x_ticks.clear();
		for (int j = 0; j < wEWM[0].size(); j++) {//ָ���������0<j<=12(�öδ��벻������if(0==i)��)
			x_t.push_back(EMS[0][j]);//ѡȡ��һ��ָ���������Ϊ�ο�
		}
		if (0 == i) {//��һ�ν������򣬺�����ͼ���մ����н��л�ͼ�������ڶ�ȡ������Ϊͬһ���͸�ʽ��
			//ָ�����EMS����
			visited = vector<int>(x_t.size(), 0);
			sortList = vector<int>(x_t.size(), 0);
			x_sort = vector<double>(x_t.size(), 0.0);
			for (int j = 0; j < x_t.size(); j++) {
				int it = 0;
				while (visited[it]) {
					if (it >= x_t.size())
						it = 0;
					else
						it++;
				}
				x_sort[j] = x_t[it];
				for (int k = 0; k < x_t.size(); k++) {
					if (x_sort[j] >= x_t[k] && visited[k] != 1) {//��������
						x_sort[j] = x_t[k];
						sortList[j] = k;
					}
				}
				visited[sortList[j]] = 1;
			}
		}
		//x_sort = x_t;
		for (int j = 0; j < wEWM[0].size(); j++) {//ָ���������0<j<=12
			x_label.push_back(j);
			wEWM_t.push_back(wEWM[i][sortList[j]]);
			weight_t(j) = wEWM[i][sortList[j]];
			xticks_t = to_string(EMS[0][sortList[j]]);
			x_ticks.push_back(IndexName[0][sortList[j]][0] + ":" + xticks_t.substr(0, 7-IndexName[0][sortList[j]][0].size())); //��֤string������ȣ��������ᱨ��(�öδ��벻������if(0 == i)��)
			
		}
		//���ݷ���
		meanValue = weight_t.mean();//ƽ��ֵ
		normValue = (weight_t - meanValue * VectorXd::Ones(weight_t.size())).norm();//��׼��
		minValue = weight_t.minCoeff();//��Сֵ
		maxValue = weight_t.maxCoeff();//���ֵ

		dataLabel.append(to_string(dataLength[i][0]) + ":" + to_string(meanValue) + "," + to_string(normValue));

		if (plotCount10++ >= 10)
			plotType[0] = 'o';
		plt::named_plot(dataLabel.c_str(), x_label, wEWM_t, plotType);

		dataLabel.clear();
		wEWM_t.clear();
		x_t.clear();
	}
	//vector<string> s;
	//for (int si = 0; si < x_label.size(); si++)
	//	s.push_back(to_string(x_label[si] + 1));
	plt::xticks(x_label,x_ticks);
	plt::title(title);
	plt::xlabel("Index & Squared Root Mean - sorted");
	plt::ylabel("Weight");
	plt::legend();
	plt::grid(1);
	plt::show();

	return;
}



vector<vector<vector<double>>> dataConvert(vector<vector<double>> indicesData,string subRateM) {
	// �ȼ�    1   2   3   4   5   6   7   8   9
	// ��Ӧֵ  1/9 1/7 1/5 1/3 1   3   5   7   9
	vector<vector<double>> dataConvert_t;
	vector<vector<vector<double>>> dataConvert_r;
	vector<double> convertTabl = { 1 / 9.0,1 / 7.0,1 / 5.0,1 / 3.0,1.0,3.0,5.0,7.0,9.0 };
	vector<vector<double>> convertTabl2 = { {0,0,0.25}, {0,0.25,0.5}, {0.25,0.5,0.75}, {0.5,0.75,1}, {0.75,1,1} };
	double tmp;
	for (int i = 0; i < indicesData.size(); i++) {
		dataConvert_t.resize(indicesData[i].size());
		for (int j = 0; j < indicesData[i].size(); j++) {
			if (indicesData[i][j] > 0)
				if (0 == subRateM.compare("AHP"))
					dataConvert_t[j].push_back(convertTabl[(int)indicesData[i][j] - 1]);
				else if (0 == subRateM.compare("FDM")) {
					dataConvert_t[j].push_back(convertTabl2[(int)indicesData[i][j] -1][0]);
					dataConvert_t[j].push_back(convertTabl2[(int)indicesData[i][j] -1][1]);
					dataConvert_t[j].push_back(convertTabl2[(int)indicesData[i][j] -1][2]);
				}
				else
					dataConvert_t[j].push_back(indicesData[i][j]);
		}
		if (!dataConvert_t.empty())
			dataConvert_r.push_back(dataConvert_t);
		dataConvert_t.clear();
	}
	return dataConvert_r;
}


vector<double> AHP(vector<vector<vector<double>>> indicesData,vector<vector<string>>& indexName) {
	vector<double> indicesVector;
	vector<double> weight,wt;
	vector<vector<double>> wtt;
	vector<vector<string>> indexName_t,indexname_o(indexName.size());
	vector<string> tt,ts;
	string str_t;
	vector<int> Lvn = { 0,3,9,12,15,16,19,22,23,26,27,(int)indicesData.size() };
	vector<int> MappingTable = { 1,2,3,10,10,10,10,20,20,20,11,11,11,12,12,13,13,13,14,14,14,21,21,22,22,22,231,31,31 };//����ָ��ӳ���
	vector<double> indicesLv1;
	vector<vector<double>> indicesLv2(indicesData.size());
	int count_indexname(0);

	//Ĭ��AHP����������άΪ1
	//��һά 28���������ڶ�άָ����������
	for (int j = 0; j < indicesData[0].size(); ++j) {
		for (int count = 0; count < Lvn.size() - 1; ++count) {
			for (int i = Lvn[count]; i < indicesData.size() && i < Lvn[count+1]; ++i) {
				for (int k = 0; k < indicesData[i][0].size(); ++k) {
					indicesVector.push_back(indicesData[i][j][k]);
				}
			}
			//һ��ָ��Ȩ��
			indexName_t.resize(Lvn[count+1]-Lvn[count]);
			for (int it = Lvn[count]; it < Lvn[count + 1]; ++it) {
				tt.assign(indexName[it].begin(), indexName[it].end());
				indexName_t[it- Lvn[count]].push_back(tt[0]);
				tt.clear();
			}
			wt = consistencyCheck(indicesVector, indexName_t);
			for (int it = 0; it < wt.size(); ++it) {
				weight.push_back(wt[it]);//ͬһάvectorֻ�ܵ�������push_back
				if (0 == j) {
					ts.assign(indexName_t[it].begin(), indexName_t[it].end());
					indexname_o[count_indexname++].push_back(ts[0]);
					ts.clear();
				}
					
			}
				
			
			indicesVector.clear();
			indexName_t.clear();
			wt.clear();
		}
		wtt.push_back(weight);
		weight.clear();
	}
	wt = vector<double>(wtt[0].size(), 0);
	int count = 0;
	for (int j = 0; j < wtt[0].size(); ++j) {
		for (int i = 0; i < wtt.size(); ++i) {
			wt[j] += wtt[i][j];
			if (wtt[i][j] > 0)
				count++;
		}
		if (count > 0) {
			wt[j] /= (double)count;
			count = 0;
		}
	}

	//����ָ���Ȩ
	weight.resize(wt.size(),-1);
	//weight.resize(19);
	count = 0;
	for (int it = 0; it < MappingTable.size(); ++it) {
		if (MappingTable[it] > 0 && MappingTable[it] < 10)
			indicesLv1.push_back(wt[it]);
		else if (MappingTable[it] > 0 && MappingTable[it] % 10 == 0)
			indicesLv2[MappingTable[it] / 10 - 1].push_back(wt[it]);
		else if (MappingTable[it] > 0 && MappingTable[it] / 10 < 10){
			if (indicesLv2[MappingTable[it] / 10 - 1].empty())//����ָ�겻����ʱ��Ĭ����Ϊ1
				weight[count++] = (wt[it] * indicesLv1[MappingTable[it] / 10 - 1] * 1);
			else
				weight[count++] = (wt[it] * indicesLv1[MappingTable[it] / 10 - 1] * indicesLv2[MappingTable[it] / 10 - 1][MappingTable[it] % 10 - 1]);
		}
		else {//һ����ָ���ڣ�����ָ�겻�ڣ�����λ�����
			weight[count++] = 1 * indicesLv1[MappingTable[it] / 100 - 1] * indicesLv2[MappingTable[it] / 100 - 1][MappingTable[it] / 10 % 10 - 1];
		}
	}
	indexName.clear();
	indexName = indexname_o;
	while (indexName.back().empty())
		indexName.pop_back();

	while (weight.back() <= 0)
		weight.pop_back();
	

	return weight;

}

vector<double> consistencyCheck(vector<double> indicesMatrix,vector<vector<string>>& indexName) {
	vector<vector<string>> indexName_t(indexName.size()+1);
	vector<string> tmp;
	string tt,tt2;
	int tag = 0, count = 0, n = 0;
	regex reg("��(.{4,24})��");
	smatch name_t;

	vector<int> xy(2,0);//��������
	VectorXd RI(15);
	double CI = 0, CR = 0;

	vector<double> weight_r;
	if (indicesMatrix.empty() || indexName.empty())
		return weight_r;

	//����ָ��
	for (int i = 0; i < indexName.size(); i++) {
		tmp.assign(indexName[i].begin(), indexName[i].end());
		tt = tmp[0].c_str();//vector->string
		auto pos = tt.cbegin();
		for (; regex_search(pos, tt.cend(), name_t, reg); pos = name_t.suffix().first) {
			//cout << name_t.str() << endl;
			for (int it = 0; it < indexName_t.size(); ++it) {
				tmp.assign(indexName_t[it].begin(), indexName_t[it].end());
				tt2 = tmp[0].c_str();
				if (0 == tt2.compare(name_t.str())) {
					tag = 1;
					break;
				}
			}
			if (0 == tag) {
				indexName_t[count++].push_back(name_t.str());//string->vector
			}
			tag = 0;
		}
	}
	for (int ni = 0; ni < indexName_t.size(); ++ni)
		if (!indexName_t[ni].empty())
			n++;
		else
			break;

	VectorXcd weight(n);
	MatrixXd judgeMatrix = MatrixXd::Identity(n, n);
	//���ɾ���
	for (int i = 0; i < indexName.size(); i++) {
		count = 0;
		tmp.assign(indexName[i].begin(), indexName[i].end());
		tt = tmp[0].c_str();
		auto pos = tt.cbegin();
		for (; regex_search(pos, tt.cend(), name_t, reg); pos = name_t.suffix().first) {
			for (int it = 0; it < indexName_t.size(); ++it) {
				tmp.assign(indexName_t[it].begin(), indexName_t[it].end());
				tt2 = tmp[0].c_str();
				if (0==tt2.compare(name_t.str()))
					xy[count++] = it;
			}
		}
		judgeMatrix(xy[0], xy[1]) = indicesMatrix[i];
		judgeMatrix(xy[1], xy[0]) = 1.0/indicesMatrix[i];
	}
	//cout << judgeMatrix << endl;

	//�����ֽ�
	EigenSolver<MatrixXd> es;
	es.compute(judgeMatrix, true);

	MatrixXd::Index maxrow;
	//cout << es.eigenvalues().rowwise().norm() << endl;
	complex<double> maxRoot = es.eigenvalues().rowwise().norm().maxCoeff(&maxrow);
	//cout << abs(maxRoot) << endl;
	MatrixXcd MRvectors(es.eigenvalues().rows(), 1);
	for (int mri = 0; mri < es.eigenvectors().rows(); ++mri)
		MRvectors(mri, 0) = es.eigenvectors()(mri, maxrow);//��ȡ�������ֵ����������
	RI << 0, 0, 0.52, 0.89, 1.12, 1.26, 1.36, 1.41, 1.46, 1.49, 1.52, 1.54, 1.56, 1.58, 1.59;

	if (maxRoot.real() > judgeMatrix.rows() && judgeMatrix.rows() > 2) {
		CI = (maxRoot.real() - judgeMatrix.rows()) / (judgeMatrix.rows() - 1);
		CR = CI / RI(judgeMatrix.rows() - 1);
		if (CR < 0.1)
			weight = MRvectors.array() / (MRvectors.sum() * VectorXd::Ones(judgeMatrix.rows(), 1)).array();
		else
			weight = VectorXd::Zero(judgeMatrix.rows(), 1);
	}
	else {
		weight = MRvectors.array() / (MRvectors.sum() * VectorXd::Ones(judgeMatrix.rows(), 1)).array();
	}
	//cout << weight << endl;
	// VectorXd -> vector<double>
	for (int jt = 0; jt < weight.size(); ++jt)
		weight_r.push_back(abs(weight(jt)));

	//ָ�꼰Ȩ���������
	indexName.clear();
	indexName.resize(indexName_t.size());
	count = 0;
	while (!indexName_t[count].empty() && count<indexName_t.size()) {
		tmp.assign(indexName_t[count].begin(), indexName_t[count].end());
		indexName[count].push_back(tmp[0]);
		count++;
	}
	while (indexName.back().empty()) {
		indexName.pop_back();
	}

	return weight_r;
}


vector<double> FDM(vector<vector<vector<double>>> indicesData, vector<vector<string>>& indexName) {
	MatrixXd minMatrix(indicesData.size(), indicesData[0].size());
	MatrixXd proMatrix(indicesData.size(), indicesData[0].size());
	MatrixXd maxMatrix(indicesData.size(), indicesData[0].size());
	MatrixXd x, y, z, ub, lb, Db;
	double alpha = 0.5, sigma = 0.5;
	vector<int> Lvn = { 0,3,7,10,11,14,17,20,23,25,28,29,31 };
	vector<int> MappingTable = { 1,2,3,10,10,10,10,20,20,20,30,11,11,11,12,12,12,13,13,13,14,14,14,21,21,22,22,22,231,31,31 };//����ָ��ӳ���

	vector<double> wt,weight,indicesLv1;
	vector<vector<double>> indicesLv2(indicesData.size());
	int count;
	//��һά��ָ��
	//�ڶ�ά������


	for (int count = 0; count < Lvn.size()-1; ++count) {
		//ÿ�γ�ʼ������
		minMatrix = MatrixXd::Zero(indicesData.size(), indicesData[0].size());
		proMatrix = MatrixXd::Zero(indicesData.size(), indicesData[0].size());
		maxMatrix = MatrixXd::Zero(indicesData.size(), indicesData[0].size());
		for (int i = Lvn[count]; i < indicesData.size() && i<Lvn[count+1]; ++i) {
			for (int j = 0; j < indicesData[i].size(); ++j) {
				minMatrix(i - Lvn[count], j) = indicesData[i][j][0];
				proMatrix(i - Lvn[count], j) = indicesData[i][j][1];
				maxMatrix(i - Lvn[count], j) = indicesData[i][j][2];
			}
		}
		x = minMatrix.rowwise().minCoeff();
		y = proMatrix.rowwise().prod();
		for (int yi = 0; yi < y.rows(); ++yi)
			y(yi, 0) = pow(y(yi, 0), 1.0 / proMatrix.cols());
		z = maxMatrix.rowwise().maxCoeff();
		ub = z - alpha * (z - y);
		lb = x - alpha * (y - x);
		Db = sigma * ub + (1 - sigma) * lb;
		Db = Db.array() / (Db.sum() * MatrixXd::Ones(Db.rows(), 1)).array();
		for (int i = 0; i < Db.rows(); ++i) {
			if (Db(i) > 0)
				wt.push_back(Db(i));
		}
		
	}

	//����ָ���Ȩ
	weight.resize(wt.size(), -1);
	//weight.resize(19);
	count = 0;
	for (int it = 0; it < MappingTable.size(); ++it) {
		if (MappingTable[it] > 0 && MappingTable[it] < 10)
			indicesLv1.push_back(wt[it]);
		else if (MappingTable[it] > 0 && MappingTable[it] % 10 == 0)
			indicesLv2[MappingTable[it] / 10 - 1].push_back(wt[it]);
		else if (MappingTable[it] > 0 && MappingTable[it] / 10 < 10) {
			if (indicesLv2[MappingTable[it] / 10 - 1].empty())//����ָ�겻����ʱ��Ĭ����Ϊ1
				weight[count++] = (wt[it] * indicesLv1[MappingTable[it] / 10 - 1] * 1);
			else
				weight[count++] = (wt[it] * indicesLv1[MappingTable[it] / 10 - 1] * indicesLv2[MappingTable[it] / 10 - 1][MappingTable[it] % 10 - 1]);
		}
		else {//һ����ָ���ڣ�����ָ�겻�ڣ�����λ�����
			weight[count++] = 1 * indicesLv1[MappingTable[it] / 100 - 1] * indicesLv2[MappingTable[it] / 100 - 1][MappingTable[it] / 10 % 10 - 1];
		}
	}

	while (weight.back() <= 0)
		weight.pop_back();
	while (indexName.back().empty()) {
		indexName.pop_back();
	}
	return weight;
}


bool PyCvxpyInputData(
	//��������Ȩ�غ�����
	string submethods, 
	vector<double> subweight, 
	vector<vector<double>> subindicesData, 
	//�͹�����Ȩ�غ�����
	string objmethods,
	vector<double> objweight,
	vector<vector<double>> objindicesData,
	//���ؽ��
	vector<vector<double>> &A,
	vector<double> &b) {


	MatrixXd SubMatrix, ObjMatrix, subWeightM, objWeightM, At, bt;

	if (subweight.empty() || subindicesData.empty() || objweight.empty() || objindicesData.empty()) {
		cout << "Empty Input" << endl;
		return false;
	}
	//һάָ�꣬��ά����
	if (subweight.size() != subindicesData.size() || objweight.size() != objindicesData.size() || subweight.size() != objweight.size()) {
		cout << "subjective and objective weight & indicesData is not Aligned!" << endl;
		return false;
	}
	//ѡȡ�ٵ����������о����ж���
	int dataLen = subindicesData[0].size() < objindicesData[0].size() ? subindicesData[0].size() : objindicesData[0].size();

	SubMatrix.resize(dataLen, subindicesData.size());
	ObjMatrix.resizeLike(SubMatrix);
	subWeightM.resize(subindicesData.size(), 1);
	objWeightM.resizeLike(subWeightM);
	//vector->Matrix
	for (int i = 0; i < subindicesData.size(); ++i) {
		subWeightM(i, 0) = subweight[i];
		objWeightM(i, 0) = objweight[i];
		for (int j = 0; j < dataLen; ++j) {
			SubMatrix(i, j) = subindicesData[i][j];
			ObjMatrix(i, j) = objindicesData[i][j];
		}
	}

	subWeightM = (subWeightM.array() / (MatrixXd::Ones(objWeightM.rows(), 1) * subWeightM.colwise().sum()).array()).matrix() ;
	At = SubMatrix + ObjMatrix;
	bt = (SubMatrix * subWeightM + ObjMatrix * objWeightM) * MatrixXd::Ones(SubMatrix.cols(), 1);
	
	A.resize(dataLen, vector<double>(bt.rows(),0));
	for (int it = 0; it < bt.rows(); ++it) {//ָ��
		b.push_back(bt(it, 0));
		for (int jt = 0; jt < dataLen; ++jt) {//����
			A[jt][it] = At(it, jt);
		}
	}

	return true;
}