#ifndef _parameter_H_
#define _parameter_H_
#include <vector>
using namespace std;

#ifndef PI
#define PI 3.1415926
#endif
class parameters {
public:
    long long fc;
    long long Bw;//bandwidth
    long long Fs;//sample frequency;
    long long speedOfLight;
    int phi;
    double PRI;
    int Nr;
    int Na;
    int N_tx;
    double dd_r;

public:
    parameters();
    virtual vector<double> windowingNr() = 0;//
    virtual vector<double> windowingNa() = 0;//
	float dr();
	float dv();
};


class para1 :public parameters {
public:
    virtual vector<double> windowingNr();
    virtual vector<double> windowingNa();
};

#endif
