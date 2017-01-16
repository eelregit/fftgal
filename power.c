#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "power.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


static double pow2(double x)
{
    return x*x;
}
static double pow3(double x)
{
    return x*x*x;
}

double amp(double x[3])
{
    return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

double *hat(double los[3])
{
    static double loshat[3];
    double losamp = amp(los);
    if(losamp > 1e-7){
        loshat[0] = los[0] / losamp;
        loshat[1] = los[1] / losamp;
        loshat[2] = los[2] / losamp;
    }
    else{
        loshat[0] = 0.;
        loshat[1] = 0.;
        loshat[2] = 0.;
    }
    return loshat;
}


void rsd(double *x, double *y, double *z, double *vx, double *vy, double *vz,
        long long int Np3, double **xd, double **yd, double **zd,
        double los[3], double aH)
{
    *xd = (double *)malloc(Np3 * sizeof(double)); assert(*xd!=NULL);
    *yd = (double *)malloc(Np3 * sizeof(double)); assert(*yd!=NULL);
    *zd = (double *)malloc(Np3 * sizeof(double)); assert(*zd!=NULL);
    double *loshat = hat(los);
    double aHinv = 1 / aH;
    clock_t t = clock();
    for(long long int p=0; p<Np3; ++p){
        double vdotlos = loshat[0]*vx[p] + loshat[1]*vy[p] + loshat[2]*vz[p];
        (*xd)[p] = x[p] + aHinv * vdotlos * loshat[0];
        (*yd)[p] = y[p] + aHinv * vdotlos * loshat[1];
        (*zd)[p] = z[p] + aHinv * vdotlos * loshat[2];
    }
    fprintf(stderr, "rsd() %.3fs to kick %lld particles along los[]={%.2f,%.2f,%.2f}\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3, los[0], los[1], los[2]);
}


int Pl(fftgal_t *fg, double dK, double los[3], char *output)
{
    int Ng = fg->Ng;
    double *loshat = hat(los);
    double dKinv = 1 / dK;
    double KF = 2 * M_PI / fg->L;
    int Nb = (int)floor(sqrt(3)*(Ng/2) * KF * dKinv) + 1;
    double *K = (double *)malloc(sizeof(double) * Nb); assert(K!=NULL);
    double *P0 = (double *)malloc(sizeof(double) * Nb); assert(P0!=NULL);
    double *P2 = (double *)malloc(sizeof(double) * Nb); assert(P2!=NULL);
    double *P4 = (double *)malloc(sizeof(double) * Nb); assert(P4!=NULL);
    double *P6 = (double *)malloc(sizeof(double) * Nb); assert(P6!=NULL);
    long int *N = (long int *)malloc(sizeof(long int) * Nb); assert(N!=NULL);
    for(int b=0; b<Nb; ++b){
        K[b] = P0[b] = P2[b] = P4[b] = P6[b] = 0.;
        N[b] = 0;
    }
    clock_t t = clock();
    for(int i=0; i<Ng; ++i){
        double Kvec[3];
        Kvec[0] = KF * remainder(i, Ng);
        for(int j=0; j<Ng; ++j){
            Kvec[1] = KF * remainder(j, Ng);
            for(int k=(i==0 && j==0); k<=Ng/2; ++k){ /* skip Kvec[]={0,0,0} */
                Kvec[2] = KF * k;
                double Kamp = amp(Kvec);
                int b = (int)floor(Kamp * dKinv);
                int count = 1 + (2*k%Ng > 0);
                double delta2 = pow2(F(fg,i,j,2*k)) + pow2(F(fg,i,j,2*k+1));
                delta2 *= count;
                double mu2 = pow2((Kvec[0]*loshat[0] + Kvec[1]*loshat[1] + Kvec[2]*loshat[2])
                        / Kamp);
                K[b] += count * Kamp;
                P0[b] += delta2;
                P2[b] += delta2 * (1.5*mu2 - 0.5);
                P4[b] += delta2 * ((4.375*mu2 - 3.75)*mu2 + 0.375);
                P6[b] += delta2 * (((14.4375*mu2 - 19.6875)*mu2 + 6.5625)*mu2 - 0.3125);
                N[b] += count;
            }
        }
    }
    fprintf(stderr, "Pl() %.3fs to bin a %d^3 grid\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Ng);

    long int Ntot = 0;
    for(int b=0; b<Nb; ++b){
        double V = pow3(fg->L);
        K[b] *= 1. / N[b];
        P0[b] *= 1. / (V * N[b]);
        P2[b] *= 5. / (V * N[b]);
        P4[b] *= 9. / (V * N[b]);
        P6[b] *= 13. / (V * N[b]);
        Ntot += N[b];
    }
    assert(Ntot == Ng*Ng*Ng-1);

    FILE *fp = fopen(output, "w"); assert(fp!=NULL);
    fprintf(fp, "# Nb Ng Np3 L dK los[3]\n");
    fprintf(fp, "%d %d %lld %f %f %f %f %f\n",
            Nb, Ng, fg->Np3, fg->L, dK, los[0], los[1], los[2]);
    fprintf(fp, "# K P0 P2 P4 P6 N\n");
    for(int b=0; b<Nb; ++b)
        fprintf(fp, "%e %e % e % e % e %ld\n", K[b], P0[b], P2[b], P4[b], P6[b], N[b]);
    fclose(fp);
    fprintf(stderr, "Pl() wrote to %s\n", output);

    free(K);
    free(P0);
    free(P2);
    free(P4);
    free(P6);
    free(N);

    return Nb;
}
