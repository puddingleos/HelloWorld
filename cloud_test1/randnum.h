#pragma once
#ifndef _RANDNUM_H_
#define _RANDNUM_H_

#include <vector>

using namespace std;
#define SRAND 2020


double GaussRand();

template<typename T>
vector<typename T> indicesDataGeneration(T refeData, T MSE, int len) {
#ifdef SRAND
	srand(SRAND);
#endif
	vector<T> indicesData;
	for (int i = 0; i < len; i++) {
		indicesData.push_back(refeData + MSE * (T)GaussRand());
	}
	return indicesData;
}


#endif
