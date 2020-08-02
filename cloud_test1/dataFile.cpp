#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "randnum.h"
#include "dataFile.h"

using namespace std;

bool dataWrite() {
	// verification to write
	string tag;
	cout << "Please Enter (yes) to begin the write:";
	cin >> tag;
	if (tag.compare("yes") != 0)
		return false;

	string Path = "";
	string filename = "indicesData";
	string filetype = ".txt";
	ofstream findicesData;
	cout << "Filename(\"indicesData\" as default): ";
	cin >> filename;

	string fid = Path.append(filename.append(filetype));
	findicesData.open(fid, ios::out | ios::trunc); // create if not exist and clear all if exist

	int indexNum = 0;
	string indexName;
	double EMS = 0;
	double refeData = 0;
	int dataLength = 0;
	vector<double> indicesData;

	cout << "Please enter the number of indicators: ";
	cin >> indexNum;
	while (indexNum--) {
		cout << "There're " << indexNum + 1 << " indices to enter.\n";
		cout << "Please enter the name of indices:";
		cin >> indexName;
		cout << "Please enter the reference, EMS and number of data:";
		cin >> refeData >> EMS >> dataLength;

		indicesData = indicesDataGeneration(refeData, EMS, dataLength);

		findicesData << indexName << "(参考值" << refeData << ",均方根" << EMS << ",数据数量" << dataLength << "): ";
		for (vector<double>::iterator it = indicesData.begin(); it != indicesData.end(); it++) {
			findicesData << *it << " ";
		}
		findicesData << "\n";

		cout << endl;
	}

	findicesData.close();

	return true;
}


bool dataRead(vector<vector<string>>& indexName, vector<vector<double>>& indicesData, vector<double>& refeData, vector<double>& EMS) {
	string Path = "";
	string filename = "indicesData";
	string filetype = ".txt";
	ifstream findicesData;
	cout << "Filename To Read(\"indicesData\" as default): ";
	cin >> filename;

	string strline;

	findicesData.open(Path.append(filename.append(filetype)),ios::in);
	if (!findicesData.is_open()) {
		cout << "file [" << filename << "] is not exist." << endl;
		return false;
	}
	

	char buff[BUFFSIZE];
	int indexNum = 0;
	while (!findicesData.eof()) {
		findicesData.getline(buff, BUFFSIZE);
		indexNum++;
	}
	findicesData.clear();//clear tag before rewind
	findicesData.seekg(0, ios::beg);

	//vector<string> indexName(indexNum);
	//vector<vector<double>> indicesData(indexNum);
	//vector<double> refeData, EMS;

	indexName.resize(indexNum-1);
	indicesData.resize(indexNum-1);
	vector<int> dataLength;
	string tmp;
	string number;
	
	int k = 0;
	while (getline(findicesData,tmp)) {
		int pos1 = tmp.find("(");
		int pos2 = tmp.find("参考值");
		int pos3 = tmp.find("均方根");
		int pos4 = tmp.find("数据数量");
		int pos5 = tmp.find(")");
		for (int i = 0; i < tmp.size(); i++) {
			if (pos1!=tmp.npos && i < pos1) 
				number.push_back(tmp[i]);
			else if(pos1 != tmp.npos && i == pos1) {
				indexName[k].push_back(number);
				number.clear();
			}
			else if (pos2 != tmp.npos && pos3 != tmp.npos && i>pos2 && i < pos3) {
				if (tmp[i] > '0'-1 && tmp[i] < '9'+1 || tmp[i] == '-' || tmp[i] == '.' || tmp[i] == 'e' || tmp[i] == 'E')
					number.push_back(tmp[i]);
			}
			else if (pos3 != tmp.npos && i == pos3) {
				refeData.push_back(atof(number.c_str()));
				number.clear();
			}
			else if (pos4 != tmp.npos && i < pos4) {
				if (tmp[i] > '0'-1 && tmp[i] < '9'+1 || tmp[i] == '-' || tmp[i] == '.' || tmp[i] == 'e' || tmp[i] == 'E')
					number.push_back(tmp[i]);
			}
			else if (pos4 != tmp.npos && i == pos4) {
				EMS.push_back(atof(number.c_str()));
				number.clear();
			}
			else if (pos5 != tmp.npos && i < pos5) {
				if (tmp[i] > '0'-1 && tmp[i] < '9'+1 || tmp[i] == '-' || tmp[i] == '.' || tmp[i] == 'e' || tmp[i] == 'E')
					number.push_back(tmp[i]);
			}
			else if (pos5 != tmp.npos && i==pos5) {
				dataLength.push_back(stoi(number.c_str()));
				number.clear();
			}
			else {
				if (tmp[i] > '0'-1 && tmp[i] < '9'+1 || tmp[i] == '-' || tmp[i] == '.' || tmp[i] == 'e' || tmp[i] == 'E')
					number.push_back(tmp[i]);
				if (tmp[i] == ' ' && !number.empty()) {
					indicesData[k].push_back(atof(number.c_str()));
					number.clear();
				}
			}
		}
		k++;
	}

	findicesData.close();
	return true;
}