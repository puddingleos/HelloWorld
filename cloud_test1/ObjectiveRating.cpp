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


//标准化数据
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
			indicesData[i][j] = abs((abs(indicesData[i][j] - refeData[i])-lowBound) / (upBound - lowBound) - dataType[i]);
		}
	}
	return true;
}

//熵权法赋权EWM
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

//灰度关联法GRA
vector<double> GreyRelationalAnalysis(vector<vector<double>> indicesData) {
	vector <double> weight(indicesData.size(), 0);
	vector<double> maxValue, minValue;
	for (vector<double>::size_type it = 0; it < indicesData.size(); it++) {
		for (vector<double>::size_type jt = 0; jt < indicesData[it].size(); jt++) {
			indicesData[it][jt] -= indicesData[it][indicesData[it].size()-1];// 参考值修正评价矩阵（参考值一般选第一列或最后一列）
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

// 主成分分析法
vector<double> PrincipleComponentAnalysis(vector<vector<double>> indicesData) {
	vector<double> weight(indicesData.size(), 0);
	vector<int> dataLength;
	vector<vector<double>> indicesData_t;
	// 以最小数据数量截取矩阵
	for (vector<double>::size_type it = 0; it < indicesData.size(); it++)
		dataLength.push_back(indicesData[it].size());
	int minLength = dataLength[min_element(dataLength.begin(), dataLength.end()) - dataLength.begin()];

	//调用Eigen库函数进行矩阵运算
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
	//2. 协方差
	MatrixXd CovMatrix;
	CovMatrix = indicesMatrix * indicesMatrix.transpose();
#ifdef DEBUG
	cout << CovMatrix << endl;
#endif
	//3. 特征分解
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
	//4. 特征根筛选
	int i = 0;
	double sum_t = 0;
	while (sum_t / eigvalueSortDescend.sum() < 0.85 && i < eigvalueSortDescend.size())
		sum_t += eigvalueSortDescend(i++);
	//5. 生成新空间
	MatrixXcd eigVectorNew = MatrixXcd::Zero(es.eigenvalues().size(), i);
	MatrixXd indicesDataSelected = MatrixXd::Zero(es.eigenvalues().size(), i);
	for (int j = 0; j < i; j++) {
		eigVectorNew.col(j) << es.eigenvectors().col(sortList(j));
		indicesDataSelected.col(j) << CovMatrix.col(sortList(j));
	}
	//cout << "eigVectorNew: " << eigVectorNew.transpose() << endl;
	//cout << "indicesDataSelected: " << indicesDataSelected << endl;
	MatrixXd T = (eigVectorNew.transpose() * indicesDataSelected).cwiseAbs2().real();

	////6. 计算指标权重
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

	// 以最小数据数量截取矩阵
	for (vector<double>::size_type it = 0; it < normIndicesData.size(); it++)
		dataLength.push_back(normIndicesData[it].size());
	int minLength = dataLength[min_element(dataLength.begin(), dataLength.end()) - dataLength.begin()];
	//调用Eigen库函数进行矩阵运算
	MatrixXd indicesMatrix((int)dataLength.size(), minLength);
	for (vector<double>::size_type it = 0; it < normIndicesData.size(); it++) {
		weight_t(it, 0) = weight[it];
		for (vector<double>::size_type jt = 0; jt < minLength; jt++) {
			indicesMatrix(it, jt) = normIndicesData[it][jt];
		}
	}
	//1. 标准矩阵赋权
	indicesMatrix = indicesMatrix.array() * (weight_t * MatrixXd::Ones(1, minLength)).array();
	//2. 正负理想解距离
	MatrixXd PostiveIdeaSolution = indicesMatrix.rowwise().maxCoeff();
	MatrixXd NegativeIdeaSolution = indicesMatrix.rowwise().minCoeff();
	MatrixXd Dp = (indicesMatrix - PostiveIdeaSolution * MatrixXd::Ones(1, minLength)).colwise().norm();//求均方根之和
	MatrixXd Dn = (indicesMatrix - NegativeIdeaSolution * MatrixXd::Ones(1, minLength)).colwise().norm();
	//3. 正理想解的贴近度(负理想解=1-正理想解)
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
	VectorXd weight_t((int)wEWM.size());//调用Eigen库矩阵运算
	double meanValue(0.0), normValue(0.0), minValue(0.0), maxValue(0.0);

	int plotCount10 = 0;
	string plotType = "*--";
	for (int j = 0; j < wEWM[0].size(); j++) {//指标数据序号0<j<=12
		for (int i = 0; i < wEWM.size(); i++) {//文件数量（场景数量）序号 0<i<10
			wEWM_t.push_back(wEWM[i][j]);
			weight_t(i) = wEWM[i][j];
			x_t.push_back((double)dataLength[i][0]);//选取第一个指标的数据作为参考
		}
		for (int k = 0; k < IndexName[0][j].size(); k++)
			IndexName_t.append(IndexName[0][j][k]);

		//数据分析
		meanValue = weight_t.mean();//平均值
		normValue = (weight_t - meanValue * VectorXd::Ones(weight_t.size())).norm();//标准差
		minValue = weight_t.minCoeff();//最小值
		maxValue = weight_t.maxCoeff();//最大值

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
	VectorXd weight_t((int)wEWM.size());//调用Eigen库矩阵运算
	double meanValue(0.0), normValue(0.0), minValue(0.0), maxValue(0.0);

	int plotCount10 = 0;
	string plotType = "*--";

		
	for (int i = 0; i < wEWM.size(); i++) {//文件数量（场景数量）序号 0<i<10
		x_label.clear();
		x_ticks.clear();
		for (int j = 0; j < wEWM[0].size(); j++) {//指标数据序号0<j<=12(该段代码不可移入if(0==i)内)
			x_t.push_back(EMS[0][j]);//选取第一个指标的数据作为参考
		}
		if (0 == i) {//第一次进行排序，后续画图依照此序列进行画图（适用于读取的数据为同一类型格式）
			//指标根据EMS排序
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
					if (x_sort[j] >= x_t[k] && visited[k] != 1) {//升序排序
						x_sort[j] = x_t[k];
						sortList[j] = k;
					}
				}
				visited[sortList[j]] = 1;
			}
		}
		//x_sort = x_t;
		for (int j = 0; j < wEWM[0].size(); j++) {//指标数据序号0<j<=12
			x_label.push_back(j);
			wEWM_t.push_back(wEWM[i][sortList[j]]);
			weight_t(j) = wEWM[i][sortList[j]];
			xticks_t = to_string(EMS[0][sortList[j]]);
			x_ticks.push_back(IndexName[0][sortList[j]][0] + ":" + xticks_t.substr(0, 7-IndexName[0][sortList[j]][0].size())); //保证string长度相等，否则后面会报错(该段代码不可移入if(0 == i)内)
			
		}
		//数据分析
		meanValue = weight_t.mean();//平均值
		normValue = (weight_t - meanValue * VectorXd::Ones(weight_t.size())).norm();//标准差
		minValue = weight_t.minCoeff();//最小值
		maxValue = weight_t.maxCoeff();//最大值

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
	// 等级    1   2   3   4   5   6   7   8   9
	// 对应值  1/9 1/7 1/5 1/3 1   3   5   7   9
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
					dataConvert_t[j][0] = convertTabl2[(int)indicesData[i][j] -1][0];
					dataConvert_t[j][1] = convertTabl2[(int)indicesData[i][j] -1][1];
					dataConvert_t[j][2] = convertTabl2[(int)indicesData[i][j] -1][2];
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


vector<double> AHP(vector<vector<vector<double>>> indicesData,vector<vector<string>> indexName) {
	vector<double> indicesVector;
	vector<double> weight(3,0);
	vector<vector<string>> indexName_t;
	vector<string> tt;
	stringstream ss;
	string str_t;
	int Lv1 = 3, Lv2 = 8, Lv3 = 19;
	//默认AHP输入矩阵第三维为1
	for (int i = 0; i < indicesData.size(); ++i) {
		for (int j = 0; j < indicesData[i].size() && j< Lv1; ++j) {
			for (int k = 0; k < indicesData[i][j].size(); ++k) {
				indicesVector.push_back(indicesData[i][j][k]);
			}
		}
		//一级指标权重
		indexName_t.resize(Lv1);
		for (int it = 0; it < Lv1; ++it) {
			tt.assign(indexName[it].begin(), indexName[it].end());
			indexName_t[it].push_back(tt[0]);
			tt.clear();
		}
		weight = consistencyCheck(indicesVector, indexName_t);
		indicesVector.clear();
		indexName_t.clear();
	}

	return weight;

}

vector<double> consistencyCheck(vector<double> indicesMatrix,vector<vector<string>> indexName) {
	VectorXcd weight(indicesMatrix.size());
	MatrixXd judgeMatrix = MatrixXd::Identity(indicesMatrix.size(), indicesMatrix.size());
	vector<vector<string>> indexName_t(indexName.size());
	vector<string> tmp;
	string tt,tt2;
	int tag = 0, count = 0;
	regex reg("“(.{4,24})”");
	smatch name_t;

	vector<int> xy(2,0);//矩阵坐标
	VectorXd RI(15);
	double CI = 0, CR = 0;

	vector<double> weight_r;

	//保存指标
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
	//生成矩阵
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
	cout << judgeMatrix << endl;

	//特征分解
	EigenSolver<MatrixXd> es;
	es.compute(judgeMatrix, true);

	MatrixXd::Index maxrow;
	cout << es.eigenvalues().rowwise().norm() << endl;
	complex<double> maxRoot = es.eigenvalues().rowwise().norm().maxCoeff(&maxrow);
	cout << abs(maxRoot) << endl;
	MatrixXcd MRvectors(es.eigenvalues().rows(), 1);
	for (int mri = 0; mri < es.eigenvectors().rows(); ++mri)
		MRvectors(mri, 0) = es.eigenvectors()(mri, maxrow);//获取最大特征值的特征向量
	RI << 0, 0, 0.52, 0.89, 1.12, 1.26, 1.36, 1.41, 1.46, 1.49, 1.52, 1.54, 1.56, 1.58, 1.59;

	if (maxRoot.real() > judgeMatrix.rows() && judgeMatrix.rows() > 2) {
		cout << judgeMatrix.rows() << endl;
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

	// VectorXd -> vector<double>
	for (int jt = 0; jt < weight.size(); ++jt)
		weight_r.push_back(abs(weight(jt)));

	return weight_r;
}