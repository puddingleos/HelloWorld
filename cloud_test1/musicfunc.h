#ifndef _musicfunc_H_
#define _musicfunc_H_
#include <vector>

using namespace std;

void sigGeneration(long double R0, long double v, vector<vector<long double>>& sigSpace);
void fft3(vector<float>& blackmanNr, vector<float>& blackmanNa, vector<vector<vector<long double>>>& sigSpace);
void kfft(vector<double>& pr, vector<double>& pi, int n, int k, vector<double>& fr, vector<double>& fi);

#endif
