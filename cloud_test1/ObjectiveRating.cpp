#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <numeric>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

#include "matplotlibcpp.h"
#include "ObjectiveRating.h"

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