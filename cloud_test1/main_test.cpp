#include <iostream>
#include <string>
#include <vector>
#include <complex>
//#include <opencv2/core.hpp>
//#include <opencv2/imgcodecs.hpp>
//#include <opencv2/highgui.hpp>
#include "randnum.h"
#include "matplotlibcpp.h"
#include "parameter.h"
#include "matprocess.h"
#include "dataFile.h"
#include "ObjectiveRating.h"

using namespace std;
//using namespace cv;
namespace plt = matplotlibcpp;


//#define TESTOPENCV
//#define MATPROCESS
#define TESTMATPLOTLIB
//#define SIGPROCESS

int main(int argc, char** argv) {

#ifdef TESTOPENCV
	string Path = "C:\\Users\\lgd\\Pictures\\Camera Roll\\";
	string filename = "bili_img_86334786.jpg";
	Mat image;
	image = imread(Path.append(filename), IMREAD_COLOR); // Read the file
	if (image.empty()) // Check for invalid input
	{
		cout << "Could not open or find the image" << std::endl;
		return -1;
	}
	namedWindow("Display window", WINDOW_AUTOSIZE); // Create a window for display.
	imshow("Display window", image); // Show our image inside it.
	waitKey(0); // Wait for a keystroke in the window
#endif // TESTOPENCV

#ifdef MATPROCESS
	// 写入文件
	//dataWrite();

	//读取文件
	vector<vector<string>> IndexName;
	vector<vector<double>> IndicesData;
	vector<double> refeData;
	vector<double> EMS;


	if (!dataRead(IndexName, IndicesData, refeData, EMS)) {
		cout << "Please check filename" << endl;
		return 0;
	}

	//数据标准化
	if (!dataMatrixNormalized(IndicesData, refeData, EMS, vector<int>(refeData.size(), 1))) {
		cout << "Matrix Normalized failed." << endl;
		return 0;
	}

	//客观评价方法
	vector<double> wEWM = EntropyWeightMethod(IndicesData);


	//plot hist
	plt::bar(wEWM);
	plt::show();

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
	//		for (int j = 0;k<sigSpace[k][i].size();j++)
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
				absRdmComb[i][j] += absRdm[i][j][k];
		}
		x.push_back(x_row);
		y.push_back(y_row);
	}

	plt::plot_surface(x,y,absRdmComb);
	plt::show();

#endif

#ifdef TESTMATPLOTLIB	
	std::vector<int> test_data;
	for (int i = 0; i < 20; i++) {
		test_data.push_back(i);
	}

	plt::bar(test_data);
	plt::show();
#endif
	return 0; // free all vectors after run this code

}