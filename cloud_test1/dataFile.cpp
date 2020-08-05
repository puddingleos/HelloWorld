#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "randnum.h"
#include "dataFile.h"

using namespace std;

bool dataWrite() {
	string Path = "";
	string filename;
	string filetype = ".txt";
	string configFiletype = ".cfg";
	ofstream findicesData;

	cout << "Please enter the fileName(indicesData): ";
	cin >> filename;

	//参数文件
	string fid = "";
	fid.append(Path);
	fid.append(filename);
	fid.append(filetype);
	findicesData.open(fid, ios::out | ios::trunc); // create if not exist and clear all if exist
	
	//配置文件（保存本次参数）
	string fid2 = "";
	fid2.append(Path);
	fid2.append(filename);
	fid2.append(configFiletype);
	

	
	int indexNum = 0;
	string indexName;
	double EMS = NULL;
	double refeData = NULL;
	int dataLength = NULL;
	vector<double> indicesData;

	// 判断是否选择上一次配置文件
	string tag;
	cout << "If you use the latest saved config (yes/no): ";
	cin >> tag;
	while (tag.compare("yes") && tag.compare("no")) {
		tag.clear();
		cout << "Please enter yes or no: ";
		cin >> tag;
	}
	if (!tag.compare("yes")) {
		ifstream configFile;
		configFile.open(fid2, ios::out);
		if (!configFile.is_open()) {
			cout << "config file is not exist." << endl;
			return false;
		}
		string tmp;
		while (getline(configFile, tmp)) {
			char *splitChar = strtok(const_cast<char*>(tmp.c_str()), " ");
			while (*splitChar != 0) {
				indexName.push_back(*splitChar);//指标名称
				splitChar++;
			}
			splitChar = strtok(NULL, " ");
			refeData = atof(splitChar);//参考值
			splitChar = strtok(NULL, " ");
			EMS = atof(splitChar);// 均方根
			splitChar = strtok(NULL, " ");
			dataLength = stoi(splitChar);//数据数量

			indicesData = indicesDataGeneration(refeData, EMS, dataLength);
			findicesData << indexName << "(参考值" << refeData << ",均方根" << EMS << ",数据数量" << dataLength << "): ";//写入数据
			for (vector<double>::iterator it = indicesData.begin(); it != indicesData.end(); it++) {
				findicesData << *it << " ";
			}
			indexName.clear();
			findicesData << "\n";

		}
		findicesData.close();
		configFile.close();
		cout << "files rewrite over." << endl;
	}
	else {
		ofstream configFile;
		configFile.open(fid2, ios::in | ios::trunc);//重新生成配置文件，以截断形式打开重名文件
		cout << "Please enter the number of indicators: ";
		cin >> indexNum;
		while (indexNum) {
			cout << "There're " << indexNum << " indices to enter.\n";
			cout << "Please enter the name of indices, reference, EMS and number of data:";
			cin >> indexName >> refeData >> EMS >> dataLength;

			indicesData = indicesDataGeneration(refeData, EMS, dataLength);
			configFile << indexName << " " << refeData << " " << EMS << " " << dataLength << "\n";//写入配置文件
			findicesData << indexName << "(参考值" << refeData << ",均方根" << EMS << ",数据数量" << dataLength << "): ";//写入数据
			for (vector<double>::iterator it = indicesData.begin(); it != indicesData.end(); it++) {
				findicesData << *it << " ";
			}
			findicesData << "\n";
			indexNum--;
			cout << endl;
		}

		findicesData.close();
		configFile.close();
	}
	

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