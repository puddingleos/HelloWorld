#pragma once
#ifndef _DATAFILE_H_
#define _DATAFILE_H_

#define BUFFSIZE 16000
#define VECTORSIZE 40
#ifndef PI
#define PI 3.1415926
#endif


bool csvRead(vector<vector<vector<string>>>& indexName_t, vector<vector<vector<double>>>& indicesData_t);
bool dataWrite(string filename = {});
bool dataWrite(string methods, vector<double> weight, vector<vector<string>> indexName, vector<vector<double>> indicesData);
bool dataWrite(string filename, vector<vector<double>> A, vector<double> b);
bool dataReadFiles(vector<vector<vector<string>>>& indexName_t, vector<vector<vector<double>>>& indicesData_t,
	vector<vector<double>>& refeData_t, vector<vector<double>>& EMS_t, vector<vector<int>>& dataLength_t);
bool dataReadFiles(string filename, vector<vector<vector<string>>>& indexName_t, vector<vector<double>>& weight_t, 
	vector<vector<vector<double>>>& indicesData_t);
int countFiles(string path, string filetype);

#endif // ! _DATAFILE_H_
