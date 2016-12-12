#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "power.h"


double pow2(double x)
{
    return x*x;
}
double pow3(double x)
{
    return x*x*x;
}

double amp(double x[3])
{
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

int normalize(double x[3])
{
    double xamp = amp(x);
    if(xamp > 1e-5){
        x[0] /= xamp;
        x[1] /= xamp;
        x[2] /= xamp;
        return 1;
    }
    else{
        x[0] = 0.;
        x[1] = 0.;
        x[2] = 0.;
        return 0;
    }
}


void rsd(double *x, double *y, double *z, double *vx, double *vy, double *vz,
        long long int Np3, double *xd, double *yd, double *zd,
        double los[3], double aH)
{
    xd = (double *)malloc(Np3 * sizeof(double)); assert(xd!=NULL);
    yd = (double *)malloc(Np3 * sizeof(double)); assert(yd!=NULL);
    zd = (double *)malloc(Np3 * sizeof(double)); assert(zd!=NULL);
    normalize(los);
    double aHinv = 1 / aH;
    for(long long int p=0; p<Np3; ++p){
        double vdotlos = los[1]*vx[p] + los[2]*vy[p] + los[3]*vz[p];
        xd[p] = x[p] + aHinv * vdotlos * los[1];
        yd[p] = y[p] + aHinv * vdotlos * los[2];
        zd[p] = z[p] + aHinv * vdotlos * los[3];
    }
}


int Pl(fftgal_t *fg, double dK, double los[3], double *K,
        double *P0, double *P2, double *P4, double *P6, long int *N)
{
    int Ng = fg->Ng;
    normalize(los);
    double dKinv = 1 / dK;
    double KF = 2 * M_PI / fg->L;
    int Nb = (int)floor(sqrt(3)*(Ng/2) * KF * dKinv) + 1;
    K = (double *)malloc(sizeof(double) * Nb); assert(K!=NULL);
    P0 = (double *)malloc(sizeof(double) * Nb); assert(P0!=NULL);
    P2 = (double *)malloc(sizeof(double) * Nb); assert(P2!=NULL);
    P4 = (double *)malloc(sizeof(double) * Nb); assert(P4!=NULL);
    P6 = (double *)malloc(sizeof(double) * Nb); assert(P6!=NULL);
    N = (long int *)malloc(sizeof(long int) * Nb); assert(N!=NULL);
    for(int b=0; b<Nb; ++b){
        K[b] = P0[b] = P2[b] = P4[b] = P6[b] = 0.;
        N[b] = 0;
    }
    for(int i=0; i<Ng; ++i)
    {
        double Kvec[3];
        Kvec[0] = KF * remainder(i, Ng);
        for(int j=0; j<Ng; ++j)
        {
            Kvec[1] = KF * remainder(j, Ng);
            for(int k=0; k<=Ng/2; ++k)
            {
                Kvec[2] = KF * k;
                double Kamp = amp(Kvec);
                int b = (int)floor(Kamp * dKinv);
                double delta2 = pow2(F(fg,i,j,2*k)) + pow2(F(fg,i,j,2*k+1));
                normalize(Kvec);
                double mu2 = pow2(Kvec[0]*los[0] + Kvec[1]*los[1] + Kvec[2]*los[2]);
                int count = (1 + (2*k%Ng > 0)) * (i*j*k > 0);
                K[b] += count * Kamp;
                P0[b] += count * delta2;
                P2[b] += count * delta2 * fma(1.5, mu2, 0.5);
                P4[b] += count * delta2 * fma(fma(4.375, mu2, -3.75), mu2, 0.375);
                P6[b] += count * delta2 *
                    fma(fma(fma(14.4375, mu2, -19.6875), mu2, 6.5625), mu2, -0.3125);
                N[b] += count;
            }
        }
    }

    long int Ntot = 0;
    for(int b=0; b<Nb; ++b)
    {
        double V = pow3(fg->L);
        K[b] /= N[b];
        P0[b] /= V * N[b];
        P2[b] /= V * N[b];
        P4[b] /= V * N[b];
        P6[b] /= V * N[b];
        Ntot += N[b];
    }
    assert(Ntot == Ng*Ng*Ng-1);

    return Nb;
}
