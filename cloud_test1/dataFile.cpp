#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "randnum.h"
#include "dataFile.h"
#include <io.h>

using namespace std;

bool dataWrite() {
	string Path = "";
	string filename;
	string filetype = ".txt";
	string configFiletype = ".cfg";
	ofstream findicesData;
	int N, nidx, options;
	double steps;
	string tag;
	int indexNum = 0;
	string indexName,fid,fid2;
	double EMS = NULL;
	double refeData = NULL;
	int dataLength = NULL;
	vector<double> indicesData;

	cout << "Please enter the numbers of files to generate:";
	cin >> N;
	cout << "Please enter the fileName(indicesData): ";
	cin >> filename;
	cout << "Which parameters you want to change in steps(1.refeData;2.EMS;3.dataLength): ";
	cin >> options;
	cout << "The steps for " << N << " files: ";
	cin >> steps;

	// 判断是否选择上一次配置文件
	
	cout << "If you use the latest saved config (yes/no): ";
	cin >> tag;
	while (tag.compare("yes") && tag.compare("no")) {
		tag.clear();
		cout << "Please enter yes or no: ";
		cin >> tag;
	}


	nidx = 0;
	while (nidx++<N){
		//参数文件
		fid = "";
		fid.append(Path);
		fid.append(filename);
		fid.append("_");
		fid.append(1,nidx + '0');
		fid.append(filetype);
		findicesData.open(fid, ios::out | ios::trunc); // create if not exist and clear all if exist


		indexNum = 0;
		EMS = NULL;
		refeData = NULL;
		dataLength = NULL;

		//配置文件（保存本次参数）
		fid2 = "";
		fid2.append(Path);
		fid2.append(filename);
		fid2.append(configFiletype);

		if (!tag.compare("yes")) {
			cout << "configurations:";

			ifstream configFile;
			configFile.open(fid2, ios::out);
			if (!configFile.is_open()) {
				cout << "config file is not exist." << endl;
				return false;
			}



			string tmp;
			while (getline(configFile, tmp)) {
				char* splitChar = strtok(const_cast<char*>(tmp.c_str()), " ");
				while (*splitChar != 0) {
					indexName.push_back(*splitChar);//指标名称
					splitChar++;
				}
				splitChar = strtok(NULL, " ");
				refeData = atof(splitChar) * (options == 1 ? steps * nidx : 1.0);//参考值
				splitChar = strtok(NULL, " ");
				EMS = atof(splitChar) * (options == 2 ? steps * nidx : 1.0);// 均方根
				splitChar = strtok(NULL, " ");
				dataLength = stoi(splitChar) * (options == 3 ? steps * nidx : 1.0);//数据数量

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
	}

	return true;
}



bool dataReadFiles(vector<vector<vector<string>>>& indexName_t, vector<vector<vector<double>>>& indicesData_t, 
	vector<vector<double>>& refeData_t, vector<vector<double>>& EMS_t, vector<vector<int>>& dataLength_t) {
	//string Path = "C:\\Users\\puddingleos\\source\\repos\\MatlabCpp3";
	string Path = "C:\\Users\\lgd\\source\\repos\\MatlabCpp";
	string searchFile = Path;
	string filename;
	string filetype = ".txt";
	ifstream findicesData;
	//cout << "Filename To Read(\"indicesData\" as default): ";
	//cin >> filename;

	vector<vector<string>> indexName;
	vector<vector<double>> indicesData;
	vector<double> refeData;
	vector<double> EMS;
	vector<int> dataLength;
	string tmp;
	string number;

	int ifiles = 0;
	int Nfiles = countFiles(searchFile, filetype);//文件数量

	//indexName_t.resize(Nfiles,vector<vector<string>>(0,vector<string>(0)));
	//indicesData_t.resize(Nfiles, vector<vector<double>>(0, vector<double>(0)));
	//refeData_t.resize(Nfiles, vector<double>(0));
	//EMS_t.resize(Nfiles, vector<double>(0));

	intptr_t hFile;//这个变量是intptr_t类型的，网上很多都是long类型。如果为long类型，那么在x64平台上就会出错！
	_finddatai64_t fileInfo;
	searchFile.append("\\*" + filetype);
	//searchFile.append("indicesData_1.txt");




	if ((hFile = _findfirst64(searchFile.c_str(), &fileInfo)) == -1)
		return false;
	do {
		filename.append(fileInfo.name);//获取文件名
		findicesData.open(Path + "\\" + filename, ios::in);
		filename.clear();
		if (!findicesData.is_open()) {
			cout << "file [" << filename << "] is not exist." << endl;
			return false;
		}


		char *buff;
		buff = (char*)malloc(sizeof(char) * BUFFSIZE);
		int indexNum = 0;
		while (!findicesData.eof()) {
			findicesData.getline(buff, BUFFSIZE);
			indexNum++;
		}
		findicesData.clear();//clear tag before rewind
		findicesData.seekg(0, ios::beg);

		//初始化
		indexName.clear();
		indicesData.clear();
		refeData.clear();
		EMS.clear();
		dataLength.clear();
		indexName.resize(indexNum - 1);
		indicesData.resize(indexNum - 1);


		int k = 0;
		while (getline(findicesData, tmp)) {
			int pos1 = tmp.find("(");
			int pos2 = tmp.find("参考值");
			int pos3 = tmp.find("均方根");
			int pos4 = tmp.find("数据数量");
			int pos5 = tmp.find(")");
			for (int i = 0; i < tmp.size(); i++) {
				if (pos1 != tmp.npos && i < pos1)
					number.push_back(tmp[i]);
				else if (pos1 != tmp.npos && i == pos1) {
					indexName[k].push_back(number);
					number.clear();
				}
				else if (pos2 != tmp.npos && pos3 != tmp.npos && i > pos2 && i < pos3) {
					if (tmp[i] > '0' - 1 && tmp[i] < '9' + 1 || tmp[i] == '-' || tmp[i] == '.' || tmp[i] == 'e' || tmp[i] == 'E')
						number.push_back(tmp[i]);
				}
				else if (pos3 != tmp.npos && i == pos3) {
					refeData.push_back(atof(number.c_str()));
					number.clear();
				}
				else if (pos4 != tmp.npos && i < pos4) {
					if (tmp[i] > '0' - 1 && tmp[i] < '9' + 1 || tmp[i] == '-' || tmp[i] == '.' || tmp[i] == 'e' || tmp[i] == 'E')
						number.push_back(tmp[i]);
				}
				else if (pos4 != tmp.npos && i == pos4) {
					EMS.push_back(atof(number.c_str()));
					number.clear();
				}
				else if (pos5 != tmp.npos && i < pos5) {
					if (tmp[i] > '0' - 1 && tmp[i] < '9' + 1 || tmp[i] == '-' || tmp[i] == '.' || tmp[i] == 'e' || tmp[i] == 'E')
						number.push_back(tmp[i]);
				}
				else if (pos5 != tmp.npos && i == pos5) {
					dataLength.push_back(stoi(number.c_str()));
					number.clear();
				}
				else {
					if (tmp[i] > '0' - 1 && tmp[i] < '9' + 1 || tmp[i] == '-' || tmp[i] == '.' || tmp[i] == 'e' || tmp[i] == 'E')
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

		indexName_t.push_back(indexName);
		refeData_t.push_back(refeData);
		EMS_t.push_back(EMS);
		indicesData_t.push_back(indicesData);
		dataLength_t.push_back(dataLength);
		ifiles++;

	} while (_findnext64(hFile, &fileInfo) != -1);
	return true;
}



int countFiles(string path,string filetype) {
	int count = 0;
	intptr_t hFile;//这个变量是intptr_t类型的，网上很多都是long类型。如果为long类型，那么在x64平台上就会出错！
	_finddatai64_t fileInfo;
	if (filetype[0] != '.')
		path.append("\\*." + filetype);
	else
		path.append("\\*" + filetype);
	if ((hFile = _findfirst64(path.c_str(), &fileInfo)) == -1)
		return 0;
	do {
		count++;
	} while (_findnext64(hFile, &fileInfo) != -1);
	return count;
}