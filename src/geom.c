#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "geom.h"


static double pow2(double x)
{
    return x*x;
}


void pbc(double *x, double *y, double *z, long long int Np3, double L)
{
    double Linv = 1 / L;
    clock_t t = clock();
    for(long long int p=0; p<Np3; ++p){
        x[p] -= L * floor(x[p] * Linv);
        y[p] -= L * floor(y[p] * Linv);
        z[p] -= L * floor(z[p] * Linv);
    }
    fprintf(stderr, "pbc() %.3fs to wrap around\n", (double)(clock()-t)/CLOCKS_PER_SEC);
}

long long int subbox(double *x, double *y, double *z, long long int Np3, double xyzlim[6],
        double **xsb, double **ysb, double **zsb, long long int Np3sb_max)
{
    double xmin = xyzlim[0], xmax = xyzlim[1];
    double ymin = xyzlim[2], ymax = xyzlim[3];
    double zmin = xyzlim[4], zmax = xyzlim[5];
    *xsb = (double *)malloc(Np3sb_max * sizeof(double)); assert(*xsb!=NULL);
    *ysb = (double *)malloc(Np3sb_max * sizeof(double)); assert(*ysb!=NULL);
    *zsb = (double *)malloc(Np3sb_max * sizeof(double)); assert(*zsb!=NULL);
    long long int p, Np3sb;
    clock_t t = clock();
    for(p=0, Np3sb=0; p<Np3 && Np3sb<Np3sb_max; ++p)
        if(x[p]>=xmin && x[p]<xmax && y[p]>=ymin && y[p]<ymax && z[p]>=zmin && z[p]<zmax){
            (*xsb)[Np3sb] = x[p];
            (*ysb)[Np3sb] = y[p];
            (*zsb)[Np3sb] = z[p];
            ++ Np3sb;
        }
    assert(p == Np3);
    *xsb = (double *)realloc(*xsb, Np3sb * sizeof(double)); assert(*xsb!=NULL);
    *ysb = (double *)realloc(*ysb, Np3sb * sizeof(double)); assert(*ysb!=NULL);
    *zsb = (double *)realloc(*zsb, Np3sb * sizeof(double)); assert(*zsb!=NULL);
    double percentage = 100.*Np3sb/Np3;
    fprintf(stderr, "subbox() %.3fs to pick out %lld/%lld particles (%.2f%%)\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3sb, Np3, percentage);
    if(Np3sb!=Np3 && percentage>99.)
        fprintf(stderr, "warning: subbox() thinks something is leaking\n");
    return Np3sb;
}

long long int subsphere(double *x, double *y, double *z, long long int Np3, double XYZR[4],
        double **xss, double **yss, double **zss, long long int Np3ss_max)
{
    double X = XYZR[0], Y = XYZR[1], Z = XYZR[2], R = XYZR[3];
    double R2 = R*R;
    *xss = (double *)malloc(Np3ss_max * sizeof(double)); assert(*xss!=NULL);
    *yss = (double *)malloc(Np3ss_max * sizeof(double)); assert(*yss!=NULL);
    *zss = (double *)malloc(Np3ss_max * sizeof(double)); assert(*zss!=NULL);
    long long int p, Np3ss;
    clock_t t = clock();
    for(p=0, Np3ss=0; p<Np3 && Np3ss<Np3ss_max; ++p)
        if(pow2(x[p]-X) + pow2(y[p]-Y) + pow2(z[p]-Z) < R2){
            (*xss)[Np3ss] = x[p];
            (*yss)[Np3ss] = y[p];
            (*zss)[Np3ss] = z[p];
            ++ Np3ss;
        }
    assert(p == Np3);
    *xss = (double *)realloc(*xss, Np3ss * sizeof(double)); assert(*xss!=NULL);
    *yss = (double *)realloc(*yss, Np3ss * sizeof(double)); assert(*yss!=NULL);
    *zss = (double *)realloc(*zss, Np3ss * sizeof(double)); assert(*zss!=NULL);
    double percentage = 100.*Np3ss/Np3;
    fprintf(stderr, "subsphere() %.3fs to pick out %lld/%lld particles (%.2f%%)\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3ss, Np3, percentage);
    if(Np3ss!=Np3 && percentage>99.)
        fprintf(stderr, "warning: subsphere() thinks something is leaking\n");
    return Np3ss;
}
