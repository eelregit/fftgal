#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "power.h"


double *hat(double los[3])
{
    static double loshat[3];
    double losamp = sqrt(gsl_pow_2(los[0]) + gsl_pow_2(los[1]) + gsl_pow_2(los[2]));
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
        long long Np3, double **xd, double **yd, double **zd,
        double los[3], double aH)
{
    *xd = (double *)malloc(Np3 * sizeof(double)); assert(*xd!=NULL);
    *yd = (double *)malloc(Np3 * sizeof(double)); assert(*yd!=NULL);
    *zd = (double *)malloc(Np3 * sizeof(double)); assert(*zd!=NULL);
    double *loshat = hat(los);
    double aHinv = 1 / aH;
    clock_t t = clock();
    for(long long p=0; p<Np3; ++p){
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
    double KF = 2*M_PI * fg->fold / fg->L;
    double dKinv = KF / dK;
    int Nb = (int)floor(M_SQRT3*(Ng/2) * dKinv) + 1;
    double *K = (double *)malloc(Nb * sizeof(double)); assert(K!=NULL);
    double *P0 = (double *)malloc(Nb * sizeof(double)); assert(P0!=NULL);
    double *P2 = (double *)malloc(Nb * sizeof(double)); assert(P2!=NULL);
    double *P4 = (double *)malloc(Nb * sizeof(double)); assert(P4!=NULL);
    double *P6 = (double *)malloc(Nb * sizeof(double)); assert(P6!=NULL);
    long int *N = (long int *)malloc(Nb * sizeof(long int)); assert(N!=NULL);
    for(int b=0; b<Nb; ++b){
        K[b] = P0[b] = P2[b] = P4[b] = P6[b] = 0.;
        N[b] = 0;
    }
    clock_t t = clock();
    double Kval[Ng];
    Kval[0] = 0.;
    for(int i=1; i<=Ng/2; ++i){
        Kval[i] = i;
        Kval[Ng-i] = - i;
    }
    for(int i=0; i<Ng; ++i)
    for(int j=0; j<Ng; ++j)
    for(int k=(i==0 && j==0); k<=Ng/2; ++k){ /* skip {0,0,0} */
        double Kamp = sqrt(gsl_pow_2(Kval[i]) + gsl_pow_2(Kval[j]) + gsl_pow_2(Kval[k]));
        int b = (int)floor(Kamp * dKinv);
        int count = 1 + (2*k%Ng > 0);
        double delta2 = gsl_pow_2(F_Re(fg,i,j,k)) + gsl_pow_2(F_Im(fg,i,j,k));
        delta2 *= count;
        double mu2 = gsl_pow_2((Kval[i]*loshat[0] + Kval[j]*loshat[1] + Kval[k]*loshat[2]) / Kamp);
        K[b] += count * Kamp;
        P0[b] += delta2;
        P2[b] += delta2 * (1.5*mu2 - 0.5);
        P4[b] += delta2 * ((4.375*mu2 - 3.75)*mu2 + 0.375);
        P6[b] += delta2 * (((14.4375*mu2 - 19.6875)*mu2 + 6.5625)*mu2 - 0.3125);
        N[b] += count;
    }
    fprintf(stderr, "Pl() %.3fs to bin a %d^3 grid\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Ng);

    long int Ntot = 0;
    for(int b=0; b<Nb; ++b){
        K[b] *= KF / N[b];
        P0[b] *= 1. / (fg->V * N[b]);
        P0[b] -= gsl_pow_6(fg->L) / fg->V / fg->Np3; /* HACK: subtract shot noise */
        P2[b] *= 5. / (fg->V * N[b]);
        P4[b] *= 9. / (fg->V * N[b]);
        P6[b] *= 13. / (fg->V * N[b]);
        Ntot += N[b];
    }
    assert(Ntot == Ng*Ng*Ng-1);

    FILE *fp = fopen(output, "w"); assert(fp!=NULL);
    fprintf(fp, "# Nb Ng Np3 L V fold dK los[3]\n");
    fprintf(fp, "%d %d %lld %f %g %d %f %f %f %f\n",
            Nb, Ng, fg->Np3, fg->L, fg->V, fg->fold, dK, los[0], los[1], los[2]);
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
