#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <Eigen/Eigenvalues>
//#include <pthread.h>
#include "randnum.h"
#include "matplotlibcpp.h"
//#include "parameter.h"
//#include "matprocess.h"
#include "dataFile.h"
#include "ObjectiveRating.h"
#include <regex>

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;



//#define EIGENTEST
//#define OBJCTIVE
//#define TESTMATPLOTLIB
//#define SIGPROCESS
//#define SUBJECTIVE
#define PYDATA


#define DEBUG

int main(int argc, char** argv) {

#ifdef EIGENTEST

	int N, tmp;
	vector<int> num;
	vector<int> HashTable(1000, 0);

	while (scanf("%d", &N) != EOF) {
		for (int i = 0; i < N; ++i) {
			cin >> tmp;
			HashTable[tmp] = 1;
		}

		for (int j = 0; j < 1000; ++j)
			if (HashTable[j])
				cout << j << endl;
		HashTable = vector<int>(1000, 0);
	}

#endif // TESTOPENCV

#ifdef OBJCTIVE

	//读取文件
	vector<vector<vector<string>>> IndexName;
	vector<vector<vector<double>>> IndicesData;
	vector<vector<double>> refeData, EMS;
	vector<vector<int>> dataLength;
	vector<vector<double>> rEWM, rGRA, rPCA;
	vector<vector<double>> wEWM, wGRA, wPCA;
	vector<double> x_t;

	// 写入文件
	dataWrite();

	if (!dataReadFiles(IndexName, IndicesData, refeData, EMS, dataLength)) {
		cout << "Please check filename" << endl;
		return 0;
	}
	for (int ifiles = 0; ifiles < IndexName.size(); ifiles++) {
		//数据标准化
		if (!dataMatrixNormalized(IndicesData[ifiles], refeData[ifiles], EMS[ifiles], vector<int>(refeData[ifiles].size(), 1))) {
			cout << "Matrix Normalized failed." << endl;
			return 0;
		}

		//客观评价方法
		
		wEWM.push_back(EntropyWeightMethod(IndicesData[ifiles]));
		wGRA.push_back(GreyRelationalAnalysis(IndicesData[ifiles]));
		wPCA.push_back(PrincipleComponentAnalysis(IndicesData[ifiles]));

		rEWM.push_back(EvlWithTOPSIS(wEWM[ifiles], IndicesData[ifiles]));
		rGRA.push_back(EvlWithTOPSIS(wGRA[ifiles], IndicesData[ifiles]));
		rPCA.push_back(EvlWithTOPSIS(wPCA[ifiles], IndicesData[ifiles]));

#ifdef DEBUG
		vector<double> wEWM_t = wEWM[ifiles];
		for (vector<double>::size_type ii = 0; ii < wEWM_t.size(); ii++)
			cout << wEWM_t[ii] << endl;
#endif
		//for (int kx = 0; kx < wEWM[ifiles].size(); kx++)
		//	x_t.push_back(kx + 1.0);
		////plot hist
		//plt::bar(wEWM[ifiles]);
		//plt::title("EntropyWeightMethod");
		//plt::show();

		//plt::bar(wGRA[ifiles]);
		//plt::title("GreyRelationalAnalysis");
		//plt::show();

		//plt::bar(wPCA[ifiles]);
		//plt::title("PrincipleComponentAnalysis");
		//plt::show();
		
		if (0==ifiles) {//写入txt文件
			dataWrite("EWM", wEWM[ifiles], IndexName[ifiles], IndicesData[ifiles]);
			dataWrite("GRA", wGRA[ifiles], IndexName[ifiles], IndicesData[ifiles]);
			dataWrite("PCA", wPCA[ifiles], IndexName[ifiles], IndicesData[ifiles]);
		}

		
		//plt::named_plot("EntropyWeightMethod", rEWM[ifiles],"*--");
		//plt::named_plot("GreyRelationalAnalysis", rGRA[ifiles],"*--");
		//plt::named_plot("PrincipleComponentAnalysis", rPCA[ifiles],"*--");
		//plt::title("files " + to_string(ifiles)+", "+to_string(rEWM[ifiles].size())+" datas");
		//plt::legend();
		//plt::grid(1);
		//plt::show();

		//x_t.clear();
	}

	//PlotWithDataLength(wEWM, dataLength, IndexName,"Entropy Weight Method");
	//PlotWithDataLength(wGRA, dataLength, IndexName, "Grey Relational Analysis");
	//PlotWithDataLength(wPCA, dataLength, IndexName, "Principle Component Analysis");
	
	//PlotWithEMS(wEWM, EMS, dataLength, IndexName, "Entropy Weight Method");
	//PlotWithEMS(wGRA, EMS, dataLength, IndexName, "Grey Relational Analysis");
	//PlotWithEMS(wPCA, EMS, dataLength, IndexName, "Principle Component Analysis");


#endif

#ifdef SUBJECTIVE
	vector<vector<vector<string>>> indexName_t, indexName_o;
	vector<vector<vector<double>>> indicesData_t, indicesData;
	vector<string> tmp, x_tick;
	vector<double> x_label;
	stringstream ss;
	vector<vector<string>> in_t;
	vector<double> weight1;


	csvRead(indexName_t, indicesData_t);
	//AHP
	indicesData = dataConvert(indicesData_t[0], "AHP");
	in_t = indexName_t[0];
	weight1 = AHP(indicesData, in_t);
	int i0 = in_t.size() - weight1.size();
	for (int i = i0; i < in_t.size(); ++i) {
		ss << i;
		x_tick.push_back("Cond" + ss.str());//不支持添加中文
		x_label.push_back(i);
		ss.str("");//stringsteam清空
	}

	if (dataWrite("AHP", weight1, in_t, indicesData_t[0]))
		cout << "data write complete!" << endl;

	plt::bar(x_label, weight1);
	plt::xticks(x_label, x_tick);
	plt::show();

	//FDM
	in_t = indexName_t[2];
	indicesData = dataConvert(indicesData_t[2], "FDM");
	vector<double> weight2 = FDM(indicesData, in_t);

	if (dataWrite("FDM", weight2, in_t, indicesData_t[2]))
		cout << "data write complete!" << endl;

	plt::bar(weight2);
	plt::show();


#endif

#ifdef PYDATA
	vector<vector<vector<string>>> indicesName_t;
	vector<vector<vector<double>>> indicesData_t, indicesData_r;
	vector<vector<double>> weight_t,weight_r,A;
	vector<double> b,refeSub,Range;
	vector<string> tmp,tmp2;
	string strA,strB;
	double weight_sum = 0;
	int indicesAlign = 0;
	int ifiles = 0;
	string filename = "*_data.txt";

	if (!dataReadFiles("FDM_data.txt", indicesName_t, weight_t, indicesData_t)) {
		cout << "No files exists." << endl;
	}
	if (!dataReadFiles("EWM_data.txt", indicesName_t, weight_t, indicesData_t)) {
		cout << "No files exists." << endl;
	}
	indicesAlign = indicesName_t[0].size();
	for (int i = 0; i < indicesName_t.size(); ++i) {
		if (indicesAlign > indicesName_t[i].size()) {
			indicesAlign = indicesName_t[i].size();
			ifiles = i;
		}
			
	}

	indicesData_r.resize(indicesData_t.size());
	weight_r.resize(weight_t.size());
		
	for (int j = 0; j < indicesAlign; ++j) {
		tmp.assign(indicesName_t[ifiles][j].begin(), indicesName_t[ifiles][j].end());
		strA = tmp[0].c_str();
		for (int k = 0; k < indicesName_t.size(); ++k) {
			for (int j = 0; j < indicesName_t[k].size(); ++j) {
				tmp2.assign(indicesName_t[k][j].begin(), indicesName_t[k][j].end());
				strB = tmp2[0].c_str();
				if (strA.find(strB) != string::npos || strB.find(strA) != string::npos) {
					indicesData_r[k].push_back(indicesData_t[k][j]);
					weight_r[k].push_back(weight_t[k][j]);
					//weight_sum += weight_t[k][j];
					break;
				}
			}
		}
	}
	//数据标准化
	refeSub = vector<double>(indicesData_r[0][0].size(), 6);
	Range = vector<double>(indicesData_r[0][0].size(), 5);
	if (!dataMatrixNormalized(indicesData_r[0], refeSub, Range, vector<int>(refeSub.size(), 0))) {
		cout << "Matrix Normalized failed." << endl;
		return 0;
	}
	if (!PyCvxpyInputData("FDM", weight_r[0], indicesData_r[0], "EWM", weight_r[1], indicesData_r[1], A, b)) {
		cout << "convert failed" << endl;
	}

	 filename = "C:\\Users\\lgd\\source\\repos\\MatlabCpp\\Matrix.txt";
	 if (dataWrite(filename, A, b))
		 cout << "writing successful!!" << endl;
#endif


#ifdef SIGPROCESS

	parameters* p = new para1();
 //   string path("D:/temp/data/non_static_reflector/scene1/");
 //   string filename("radar_echo_20200115201411.dat");
 //   string dataHeadFlag("ZZZZZZZZ");
 //   string dataFlag("radarEcho");
 //   unsigned long pointer = 0;
	//vector3DS_t radarcube;
	//vector3DCLD_t sigSpace(p->N_tx, vector2DCLD_t(p->Nr, vector1DCLD_t(p->Na,0)));
 //   if(!measureDataRead(path.append(filename),dataHeadFlag,dataFlag,pointer,radarcube))
	//	return 1;
	//for (int k = 0; k < sigSpace.size(); k++) {
	//	for (int i = 0; i < sigSpace[k].size(); i++) {
	//		for (int j = 0;j<sigSpace[k][i].size();j++)
	//			sigSpace[k][i][j] = radarcube[i][j][k];
	//	}
	//}
	
	

	int R0 = 12;
	int v = 20;
	int theta = 60;
	vector3DCLD_t sigSpace;
	sigGeneration((long double)R0, (long double)v, (long double)theta, sigSpace);

	//窗函数
	vector<double> blackmanNr, blackmanNa;
	blackmanNr = p->windowingNr();
	blackmanNa = p->windowingNa();
	//二维FFT
	vector3DCLD_t Rdm;
	vector3DLD_t absRdm((p->Nr) / 2, vector2DLD_t(p->Na, vector1DLD_t(p->N_tx)));
	Rdm = fft3(blackmanNr, blackmanNr, sigSpace);


	// antenna combining
	vector1DLD_t sigQuickTime;
	// abs of RD specture (Nr/2,Na,Nrx)
	for (int rdm_ri = 0; rdm_ri < (p->Nr) / 2; rdm_ri++) {
		for (int rdm_rx = 0; rdm_rx < p->N_tx; rdm_rx++) {
			for (int rdm_ra = 0; rdm_ra < p->Na; rdm_ra++) {
				absRdm[rdm_ri][rdm_ra][rdm_rx] = abs(Rdm[rdm_rx][rdm_ri][rdm_ra]);
			}
		}
	}
	vector2DLD_t absRdmComb(absRdm.size(), vector1DLD_t(absRdm[0].size(),0)),x,y;
	for (int i = 0; i <absRdm.size(); i++) {
		vector<long double>x_row, y_row;
		for (int j = 0; j < absRdm[i].size(); j++) {
			x_row.push_back(i);
			y_row.push_back(j);
			for (int k = 0; k < absRdm[i][j].size(); k++)
				absRdmComb[i][j] += log10(absRdm[i][j][k]);
		}
		x.push_back(x_row);
		y.push_back(y_row);
	}

	plt::plot_surface(x,y,absRdmComb);
	plt::show();

#endif

#ifdef TESTMATPLOTLIB	
	// Prepare data.
	int n = 5000;
	std::vector<double> x(n), y(n), z(n), w(n, 2);
	for (int i = 0; i < n; ++i) {
		x.at(i) = i * i;
		y.at(i) = sin(2 * PI * i / 360.0);
		z.at(i) = log(i);
	}

	// Set the size of output image = 1200x780 pixels
	plt::figure_size(1200, 780);

	// Plot line from given x and y data. Color is selected automatically.
	plt::plot(x, y);

	// Plot a red dashed line from given x and y data.
	plt::plot(x, w, "r--");

	// Plot a line whose name will show up as "log(x)" in the legend.
	plt::named_plot("log(x)", x, z);

	// Set x-axis to interval [0,1000000]
	plt::xlim(0, 1000 * 1000);

	// Add graph title
	plt::title("Sample figure");

	// Enable legend.
	plt::legend();
	plt::show();
	// save figure
	const char* filename = "./basic.png";
	std::cout << "Saving result to " << filename << std::endl;
	plt::save(filename);

#endif



	return 0; // free all vectors after run this code

}
