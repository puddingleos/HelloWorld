#ifndef MATPROCESS_H
#define MATPROCESS_H

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <complex>
using namespace std;

#define sample_crp 512
#define crp_frame 128
#define rx_n 4
#define PI 3.1415926

typedef vector<vector<vector<short>>> vector3DS_t;
typedef vector<vector<short>> vector2DS_t;
typedef vector<short> vector1DS_t;
typedef vector<vector<vector<long double>>> vector3DLD_t;
typedef vector<vector<long double>> vector2DLD_t;
typedef vector<long double> vector1DLD_t;
typedef vector<vector<vector<complex<long double>>>> vector3DCLD_t;
typedef vector<vector<complex<long double>>> vector2DCLD_t;
typedef vector<complex<long double>> vector1DCLD_t;



bool measureDataRead(string filename, string dataHeadFlag, string dataFlag, unsigned long pointer, vector3DS_t& radarcube);
void sigGeneration(long double R0, long double v, long double theta, vector3DCLD_t& sigSpace);
void kfft(vector1DLD_t& pr, vector1DLD_t& pi, int n, int k, vector1DLD_t& fr, vector1DLD_t& fi);
vector3DCLD_t fft3(vector<double>& blackmanNr, vector<double>& blackmanNa, vector3DCLD_t& sigSpace);

#endif // MATPROCESS_H
