#pragma once
#ifndef _DATAFILE_H_
#define _DATAFILE_H_

#define BUFFSIZE 1024

bool dataWrite();
bool dataRead(vector<vector<string>>& indexName, vector<vector<double>>& indicesData, vector<double>& refeData, vector<double>& EMS);


#endif // ! _DATAFILE_H_
