#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "geom.h"


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

double pdist2(double x1, double y1, double z1, double x2, double y2, double z2, double L)
{
    double dx = remainder(x1-x2, L);
    double dy = remainder(y1-y2, L);
    double dz = remainder(z1-z2, L);
    return dx*dx + dy*dy + dz*dz;
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

long long int subsphere(double *x, double *y, double *z, long long int Np3, double XYZRL[5],
        double **xsb, double **ysb, double **zsb, long long int Np3ss_max)
{
    double X = XYZRL[0], Y = XYZRL[1], Z = XYZRL[2], R = XYZRL[3], L=XYZRL[4];
    double R2 = R*R;
    *xsb = (double *)malloc(Np3ss_max * sizeof(double)); assert(*xsb!=NULL);
    *ysb = (double *)malloc(Np3ss_max * sizeof(double)); assert(*ysb!=NULL);
    *zsb = (double *)malloc(Np3ss_max * sizeof(double)); assert(*zsb!=NULL);
    long long int p, Np3ss;
    clock_t t = clock();
    for(p=0, Np3ss=0; p<Np3 && Np3ss<Np3ss_max; ++p)
        if(pdist2(x[p], y[p], z[p], X, Y, Z, L) < R2){
            (*xsb)[Np3ss] = x[p];
            (*ysb)[Np3ss] = y[p];
            (*zsb)[Np3ss] = z[p];
            ++ Np3ss;
        }
    assert(p == Np3);
    /* fill the remaining space of the bounding box with Np3rand = Np3sb-Np3ss,
     * with the shot noise kept the same */
    long long int Np3sb = round(6/M_PI * Np3ss);
    *xsb = (double *)realloc(*xsb, Np3sb * sizeof(double)); assert(*xsb!=NULL);
    *ysb = (double *)realloc(*ysb, Np3sb * sizeof(double)); assert(*ysb!=NULL);
    *zsb = (double *)realloc(*zsb, Np3sb * sizeof(double)); assert(*zsb!=NULL);
    static gsl_rng *rng = NULL;
    if(rng == NULL){
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_default);
    }
    for(p=Np3ss; p<Np3sb; ++p){
        double xrand, yrand, zrand;
        do{
            xrand = X + R * (2*gsl_rng_uniform(rng) - 1);
            yrand = Y + R * (2*gsl_rng_uniform(rng) - 1);
            zrand = Z + R * (2*gsl_rng_uniform(rng) - 1);
        }while(pdist2(xrand, yrand, zrand, X, Y, Z, L) < R2);
        (*xsb)[p] = xrand;
        (*ysb)[p] = yrand;
        (*zsb)[p] = zrand;
    }
    double percentage = 100.*Np3sb/Np3;
    fprintf(stderr, "subsphere() %.3fs to pick out %lld(%lld)/%lld particles (%.2f%%)\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3sb, Np3ss, Np3, percentage);
    return Np3sb;
}
