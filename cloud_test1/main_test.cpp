#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "randnum.h"
using namespace std;


int main() {

	string Path = "";
	string filename = "indicesData";
	string filetype = ".txt";
	ofstream findicesData;

	string fid = Path.append(filename.append(filetype));
	findicesData.open(fid, ios::out | ios::trunc); // create if not exist and clear all if exist
	
	string indexName;
	vector<double> indicesData;
	indexName = "距离检测结果";
    indicesData = indicesDataGeneration(5.0, 0.05, 10);
	
	findicesData << indexName << "(参考值"<< 5.0 <<"): ";
	for (vector<double>::iterator it= indicesData.begin(); it != indicesData.end(); it++) {
		findicesData << *it << " ";
	}
	findicesData << "\n";


	findicesData.close();

	return 0;
}