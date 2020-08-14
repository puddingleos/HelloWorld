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


#define EIGENTEST
//#define MATPROCESS
//#define TESTMATPLOTLIB
//#define SIGPROCESS

#define DEBUG

int main(int argc, char** argv) {

#ifdef EIGENTEST


	//MatrixXd judgeMatrix(3, 3);
	//judgeMatrix << 1, 9, 7,
	//	1.0 / 9, 1, 5,
	//	1.0 / 7, 1.0 / 5, 1;
	////cout << judgeMatrix << endl;
	////特征分解
	//EigenSolver<MatrixXd> es;
	//es.compute(judgeMatrix, true);

	//MatrixXcd::Index maxrow,maxcol;
	//complex<double> maxRoot = es.eigenvalues().rowwise().squaredNorm().maxCoeff(&maxrow);
	//MatrixXcd eigvec = es.eigenvectors();
	////std::cout << eigvec.real() << endl;
	//MatrixXcd MRvectors(eigvec.rows(),1);
	//for (int mri = 0; mri < eigvec.rows();++mri)
	//	MRvectors(mri,0) = eigvec(mri,maxrow);//获取最大特征值的特征向量


	//VectorXd RI(15);
	//double CI, CR;
	//ArrayXcd weight;
	//vector<double> weight_r;
	//RI << 0, 0, 0.52, 0.89, 1.12, 1.26, 1.36, 1.41, 1.46, 1.49, 1.52, 1.54, 1.56, 1.58, 1.59;

	//if (abs(maxRoot) > judgeMatrix.rows() && judgeMatrix.rows() > 2) {
	//	CI = (abs(maxRoot) - judgeMatrix.rows()) / (judgeMatrix.rows() - 1);
	//	CR = CI / RI(judgeMatrix.rows() - 1);
	//	//cout << abs(MRvectors.sum()) << endl;
	//	//cout << MRvectors.sum() * MatrixXd::Ones(judgeMatrix.rows(), 1) << endl;
	//	//cout << MRvectors << endl;
	//	if (CR < 0.1) {
	//		
	//		weight = MRvectors.array() / (MRvectors.sum() * MatrixXd::Ones(judgeMatrix.rows(), 1)).array();
	//	}
	//		
	//	else
	//		weight = VectorXd::Zero(judgeMatrix.rows(), 1);
	//}
	//else {
	//	weight = MRvectors.array() / (MRvectors.sum() * MatrixXd::Ones(judgeMatrix.rows(), 1)).array();
	//	
	//}

	//// VectorXd -> vector<double>
	//for (int jt = 0; jt < weight.size(); ++jt)
	//	weight_r.push_back(abs(weight(jt)));



	vector<vector<vector<string>>> indexName_t; 
	vector<vector<vector<double>>> indicesData_t,indicesData;
	csvRead(indexName_t, indicesData_t);

	indicesData = dataConvert(indicesData_t[0], "AHP");
	vector<double> weight = AHP(indicesData, indexName_t[0]);

#endif // TESTOPENCV

#ifdef MATPROCESS

	//读取文件
	vector<vector<vector<string>>> IndexName;
	vector<vector<vector<double>>> IndicesData;
	vector<vector<double>> refeData, EMS;
	vector<vector<int>> dataLength;
	vector<vector<double>> rEWM, rGRA, rPCA;
	vector<vector<double>> wEWM, wGRA, wPCA;
	vector<double> x_t;

	// 写入文件
	//dataWrite();

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

		
		plt::named_plot("EntropyWeightMethod", rEWM[ifiles],"*--");
		plt::named_plot("GreyRelationalAnalysis", rGRA[ifiles],"*--");
		plt::named_plot("PrincipleComponentAnalysis", rPCA[ifiles],"*--");
		plt::title("files " + to_string(ifiles)+", "+to_string(rEWM[ifiles].size())+" datas");
		plt::legend();
		plt::grid(1);
		plt::show();

		//x_t.clear();
	}

	//PlotWithDataLength(wEWM, dataLength, IndexName,"Entropy Weight Method");
	//PlotWithDataLength(wGRA, dataLength, IndexName, "Grey Relational Analysis");
	//PlotWithDataLength(wPCA, dataLength, IndexName, "Principle Component Analysis");
	
	//PlotWithEMS(wEWM, EMS, dataLength, IndexName, "Entropy Weight Method");
	//PlotWithEMS(wGRA, EMS, dataLength, IndexName, "Grey Relational Analysis");
	//PlotWithEMS(wPCA, EMS, dataLength, IndexName, "Principle Component Analysis");



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
