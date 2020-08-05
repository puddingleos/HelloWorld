#include <iostream>
#include <math.h>
#include <vector>
//#include "musicfunc.h"
#include "parameter.h"
#define PI 3.1415926
using namespace std;

void sigGeneration(long double R0, long double v, vector<vector<long double>>&sigSpace) {
	parameters* p = new para1;
	long double Tp = (*p).Nr / (*p).Fs;
	long double kr = (*p).Bw / Tp;
	long double dr = (*p).dr();
	long double dv = (*p).dv();
	long double c = (*p).speedOfLight;
	long double R,tr;
	vector<long double> sigQuickTime((*p).Nr, 0);
	for (int na = 0;na < (*p).Na;na++) {
		R = R0 + v * na * (*p).PRI;
		for (int nr = 0;nr < (*p).Nr;nr++) {
			tr = (nr - (*p).Nr / 2) / (*p).Fs;
			sigQuickTime[nr] = (cos(2 * PI * (tr * (kr * 2.0 * R / c - 4.0 * v * kr * R / c / c + 2.0 * v * (*p).fc / c)
				+ tr * tr * (2 * v * kr / c - 2 * v * v * kr / c / c) + ((*p).fc * (*p).fc * R / c - kr * 2 * R *
					R / c / c - (*p).phi))));
		}
		sigSpace.push_back(sigQuickTime);
		sigQuickTime.clear();
	}
}


void fft3(vector<float>& blackmanNr, vector<float>& blackmanNa, vector<vector<vector<long double>>>& sigSpace) {
	parameters* p = new para1;

}



void kfft(vector<double> &pr, vector<double> &pi, int n, int k, vector<double> &fr, vector<double> &fi)
{
	int it, m, is, i, j, nv, l0;
	double p, q, s, vr, vi, poddr, poddi;
	for (it = 0; it <= n - 1; it++)  //将pr[0]和pi[0]循环赋值给fr[]和fi[]
	{
		m = it;
		is = 0;
		for (i = 0; i <= k - 1; i++)
		{
			j = m / 2;
			is = 2 * is + (m - 2 * j);
			m = j;
		}
		fr[it] = pr[is];
		fi[it] = pi[is];
	}
	pr[0] = 1.0;
	pi[0] = 0.0;
	p = 6.283185306 / (1.0 * n);
	pr[1] = cos(p); //将w=exp(-j2pi/n)用欧拉公式表示
	pi[1] = -sin(p);

	for (i = 2; i <= n - 1; i++)  //计算pr[]
	{
		p = pr[i - 1] * pr[1];
		q = pi[i - 1] * pi[1];
		s = (pr[i - 1] + pi[i - 1]) * (pr[1] + pi[1]);
		pr[i] = p - q; pi[i] = s - p - q;
	}
	for (it = 0; it <= n - 2; it = it + 2)
	{
		vr = fr[it];
		vi = fi[it];
		fr[it] = vr + fr[it + 1];
		fi[it] = vi + fi[it + 1];
		fr[it + 1] = vr - fr[it + 1];
		fi[it + 1] = vi - fi[it + 1];
	}
	m = n / 2;
	nv = 2;
	for (l0 = k - 2; l0 >= 0; l0--) //蝶形运算
	{
		m = m / 2;
		nv = 2 * nv;
		for (it = 0; it <= (m - 1) * nv; it = it + nv)
			for (j = 0; j <= (nv / 2) - 1; j++)
			{
				p = pr[m * j] * fr[it + j + nv / 2];
				q = pi[m * j] * fi[it + j + nv / 2];
				s = pr[m * j] + pi[m * j];
				s = s * (fr[it + j + nv / 2] + fi[it + j + nv / 2]);
				poddr = p - q;
				poddi = s - p - q;
				fr[it + j + nv / 2] = fr[it + j] - poddr;
				fi[it + j + nv / 2] = fi[it + j] - poddi;
				fr[it + j] = fr[it + j] + poddr;
				fi[it + j] = fi[it + j] + poddi;
			}
	}
	for (i = 0; i <= n - 1; i++)
	{
		pr[i] = sqrt(fr[i] * fr[i] + fi[i] * fi[i]);  //幅值计算
	}
	return;
}
//————————————————
//版权声明：本文为CSDN博主「杨贵安」的原创文章，遵循CC 4.0 BY - SA版权协议，转载请附上原文出处链接及本声明。
//原文链接：https ://blog.csdn.net/yga_airspace/article/details/86688278