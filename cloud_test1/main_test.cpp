#include <iostream>
#include <string>
#include <vector>
//#include <opencv2/core.hpp>
//#include <opencv2/imgcodecs.hpp>
//#include <opencv2/highgui.hpp>
#include "randnum.h"
#include "dataFile.h"
//using namespace cv;
using namespace std;

//#define TESTOPENCV

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



	dataWrite();

	//vector<vector<string>> IndexName;
	//vector<vector<double>> IndicesData;
	//vector<double> refeData;
	//vector<double> EMS;

	//vector<vector<string>>& idn = IndexName;
	//vector<vector<double>>& idD = IndicesData;
	//vector<double>& reD = refeData;
	//vector<double>& ems = EMS;

	//if (!dataRead(idn, idD, reD, ems)) {
	//	cout << "Please check filename" << endl;
	//	return 0;
	//}
	//	
	


	return 0; // free all vectors after run this code

}