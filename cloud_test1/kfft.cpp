#include <iostream>
#include <math.h>
#include <vector>
#include <complex>
#include "matplotlibcpp.h"
#include "matprocess.h"
#include "parameter.h"


#define DEBUG
using namespace std;
namespace plt = matplotlibcpp;



void sigGeneration(long double R0, long double v, long double theta, vector3DCLD_t&sigSpace) {
    parameters* p = new para1;
    long double Tp = (long double)(*p).Nr / (long double)(*p).Fs;
    long double kr = (long double)(*p).Bw / (long double)Tp;
    long double dr = (*p).dr();
    long double dv = (*p).dv();
    long double c = (*p).speedOfLight;
    long double R,tr;
//    vector<complex<long double>> sigQuickTime((*p).Nr, 0);
    vector2DCLD_t sigspace;
    sigspace.resize((*p).Nr, vector1DCLD_t((*p).Na));

    for (int na = 0;na < (*p).Na;na++) {
        R = R0 + v * na * (*p).PRI;
        for (int nr = 0;nr < (*p).Nr;nr++) {
            tr = (nr - (*p).Nr / 2.0) / (*p).Fs;
            sigspace[nr][na] = cos(2.0 * PI * (tr * (kr * 2.0 * R / c - 4.0 * v * kr * R / c / c + 2.0 * v * (*p).fc / c)
                + tr * tr * (2.0 * v * kr / c - 2.0 * v * v * kr / c / c) + ((*p).fc * 2.0 * R / c - kr * 2.0 * R *
                    R / c / c - (*p).phi)));
        }
    }
    for (int rx_i = 0; rx_i < p->N_tx; rx_i++) {
        complex<long double> tmp = { 0,2 * PI * p->dd_r * rx_i * sin(theta * PI / 180) };
        for (int na_i = 0; na_i < p->Na; na_i++) {
            for (int nr_i = 0; nr_i < p->Nr; nr_i++) {

                sigspace[nr_i][na_i] += sigspace[nr_i][na_i] * exp(tmp);
            }
        }
        sigSpace.push_back(sigspace);
    }
}


vector3DCLD_t fft3(vector<double>& blackmanNr, vector<double>& blackmanNa, vector3DCLD_t& sigSpace) {
    parameters* p = new para1;
    vector3DCLD_t Rdm(sigSpace.size(), vector2DCLD_t((*p).Nr, vector1DCLD_t((*p).Na)));
    vector3DLD_t Rdm_r(p->N_tx, vector2DLD_t(p->Nr, vector1DLD_t(p->Na))),
            Rdm_i(Rdm_r);
    vector1DLD_t sigQuickTime_r(p->Nr,0.0), sigSlowTime_r(p->Na,0.0),sigQuickTime_i(p->Nr,0.0), sigSlowTime_i(p->Na,0.0),fr,fi;

    //range fft
    for (int rx_i = 0; rx_i < p->N_tx; rx_i++) {
        for(int na_i=0 ; na_i < p->Na; na_i++){
            for(int nr_i=0 ; nr_i < p->Nr; nr_i++){
                complex<long double> tmp = sigSpace[rx_i][nr_i][na_i]*(complex<long double>)blackmanNr[nr_i];
                sigQuickTime_r[nr_i] = (tmp.real());
                sigQuickTime_i[nr_i] = (tmp.imag());
                fr.push_back(0.0);
                fi.push_back(0.0);
            }

            kfft(sigQuickTime_r,sigQuickTime_i,p->Nr,log2(p->Nr),fr,fi);
            for (int rdm_ri = 0; rdm_ri < p->Nr; rdm_ri++) {
                Rdm_r[rx_i][rdm_ri][na_i] = sigQuickTime_r[rdm_ri];
                Rdm_i[rx_i][rdm_ri][na_i] = sigQuickTime_i[rdm_ri];
            }
            fr.clear();
            fi.clear();
        }
    }

#ifdef DEBUG
    vector2DLD_t absRdmDebug(Rdm_r[0].size(), vector1DLD_t(Rdm_r[0][0].size(), 0)), x1, y1;
    for (int i = 0; i < Rdm_r[0].size(); i++) {
        vector<long double>x_row, y_row;
        for (int j = 0; j < Rdm_r[0][i].size(); j++) {
            x_row.push_back(i);
            y_row.push_back(j);
            complex<long double> tmp = { Rdm_r[0][i][j],Rdm_i[0][i][j] };
            absRdmDebug[i][j] += log10(abs(tmp));
        }
        x1.push_back(x_row);
        y1.push_back(y_row);
    }

    plt::plot_surface(x1, y1, absRdmDebug);
    plt::show();

#endif

    // doppler fft
    for (int rx_i = 0; rx_i < p->N_tx; rx_i++) {
        for(int nr_i=0 ; nr_i < p->Nr; nr_i++){
            for(int na_i=0 ; na_i < p->Na; na_i++){
                complex<long double> tmp = {Rdm_r[rx_i][nr_i][na_i],Rdm_i[rx_i][nr_i][na_i]};
                tmp = tmp*(complex<long double>)blackmanNa[na_i];
                sigSlowTime_r[na_i] = (tmp.real());
                sigSlowTime_i[na_i] = (tmp.imag());
                fr.push_back(0.0);
                fi.push_back(0.0);
            }

            kfft(sigSlowTime_r,sigSlowTime_i,p->Na,log2(p->Na),fr,fi);
            for(int rdm_ra = 0;rdm_ra<p->Na; rdm_ra++){
                Rdm[rx_i][nr_i][rdm_ra] = { sigSlowTime_r[rdm_ra],sigSlowTime_i[rdm_ra]};
            }
            fr.clear();
            fi.clear();
        }
    }

#ifdef DEBUG
    vector2DLD_t absRdmDebug2(Rdm[0].size(), vector1DLD_t(Rdm[0][0].size(), 0)), x2, y2;
    for (int i = 0; i < Rdm[0].size(); i++) {
        vector<long double>x_row, y_row;
        for (int j = 0; j < Rdm[0][i].size(); j++) {
            x_row.push_back(i);
            y_row.push_back(j);
            absRdmDebug2[i][j] += log10(abs(Rdm[0][i][j]));
        }
        x2.push_back(x_row);
        y2.push_back(y_row);
    }

    plt::plot_surface(x2, y2, absRdmDebug2);
    plt::show();

#endif

    return Rdm;
}



void kfft(vector1DLD_t&pr, vector1DLD_t&pi, int n, int k, vector1DLD_t&fr, vector1DLD_t&fi)
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
