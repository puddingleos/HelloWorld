#include <iostream>
#include <string>
#include <vector>
#include <complex>
#include <Eigen/Eigenvalues>
#include "randnum.h"
#include "matplotlibcpp.h"
//#include "parameter.h"
//#include "matprocess.h"
#include "dataFile.h"
#include "ObjectiveRating.h"

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;


//#define EIGENTEST
//#define MATPROCESS
#define TESTMATPLOTLIB
//#define SIGPROCESS

//#define DEBUG

int main(int argc, char** argv) {



#ifdef EIGENTEST
	
	MatrixXd m(3,3);
	double tmp;
	for (int i = 0; i < m.rows(); i++) {
		for (int j = 0; j < m.cols(); j++) {
			m(i, j) = cin.get();
		}
	}
	
	cout << m << endl;
	EigenSolver<MatrixXd> es(m);
	VectorXcd eigvalue = es.eigenvalues();
	cout << eigvalue << endl;

	cout << es.eigenvectors() << endl;


#endif // TESTOPENCV

#ifdef MATPROCESS
	// 写入文件
	//dataWrite();

	//读取文件
	vector<vector<vector<string>>> IndexName;
	vector<vector<vector<double>>> IndicesData;
	vector<vector<double>> refeData;
	vector<vector<double>> EMS;


	if (!dataReadFiles(IndexName, IndicesData, refeData, EMS)) {
		cout << "Please check filename" << endl;
		return 0;
	}
	//for (int ifiles = 0; ifiles < IndexName.size(); ifiles++) {
	//	//数据标准化
	//	if (!dataMatrixNormalized(IndicesData[ifiles], refeData[ifiles], EMS[ifiles], vector<int>(refeData[ifiles].size(), 1))) {
	//		cout << "Matrix Normalized failed." << endl;
	//		return 0;
	//	}

	//	//客观评价方法
	//	vector<vector<double>> wEWM, wGRA, wPCA;
	//	wEWM.push_back(EntropyWeightMethod(IndicesData[ifiles]));
	//	wGRA.push_back(GreyRelationalAnalysis(IndicesData[ifiles]));
	//	wPCA.push_back(PrincipleComponentAnalysis(IndicesData[ifiles]));

		
#ifdef DEBUG
		vector<double> wEWM_t = wEWM[ifiles];
		for (vector<double>::size_type ii = 0; ii < wEWM_t.size(); ii++)
			cout << wEWM_t[ii] << endl;
#endif
		//plot hist
		plt::bar(IndicesData[0][0]);
		plt::title("EntropyWeightMethod");
		plt::show();

		//plt::bar(wGRA.back());
		//plt::title("GreyRelationalAnalysis");
		//plt::show();

		//plt::bar(wPCA.back());
		//plt::title("PrincipleComponentAnalysis");
		//plt::show();
	//}
	

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
	std::cout << "Saving result to " << filename << std::endl;;
	plt::save(filename);

#endif
	return 0; // free all vectors after run this code

}