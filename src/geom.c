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

long long int subbox(double *x, double *y, double *z, long long int Np3, double xyzlim[6], double bbox[6],
        double **xbb, double **ybb, double **zbb, long long int Np3sub_max)
{
    double xmin = xyzlim[0], xmax = xyzlim[1];
    double ymin = xyzlim[2], ymax = xyzlim[3];
    double zmin = xyzlim[4], zmax = xyzlim[5];
    double Vsub = (xmax-xmin) * (ymax-ymin) * (zmax-zmin);
    assert(xmin<xmax && ymin<ymax && zmin<zmax);

    double xbbmin = bbox[0], xbbmax = bbox[1];
    double ybbmin = bbox[2], ybbmax = bbox[3];
    double zbbmin = bbox[4], zbbmax = bbox[5];
    double Vbb = (xbbmax-xbbmin) * (ybbmax-ybbmin) * (zbbmax-zbbmin);
    assert(xmin>=xbbmin && xmax<=xbbmax);
    assert(ymin>=ybbmin && ymax<=ybbmax);
    assert(zmin>=zbbmin && zmax<=zbbmax);

    *xbb = (double *)malloc(Np3sub_max * sizeof(double)); assert(*xbb!=NULL);
    *ybb = (double *)malloc(Np3sub_max * sizeof(double)); assert(*ybb!=NULL);
    *zbb = (double *)malloc(Np3sub_max * sizeof(double)); assert(*zbb!=NULL);
    long long int p, Np3sub;
    clock_t t = clock();
    for(p=0, Np3sub=0; p<Np3 && Np3sub<Np3sub_max; ++p)
        if(x[p]>=xmin && x[p]<xmax && y[p]>=ymin && y[p]<ymax && z[p]>=zmin && z[p]<zmax){
            (*xbb)[Np3sub] = x[p];
            (*ybb)[Np3sub] = y[p];
            (*zbb)[Np3sub] = z[p];
            ++ Np3sub;
        }
    assert(p == Np3);

    long long int Np3bb = round(Vbb/Vsub * Np3sub);
    *xbb = (double *)realloc(*xbb, Np3bb * sizeof(double)); assert(*xbb!=NULL);
    *ybb = (double *)realloc(*ybb, Np3bb * sizeof(double)); assert(*ybb!=NULL);
    *zbb = (double *)realloc(*zbb, Np3bb * sizeof(double)); assert(*zbb!=NULL);
    static gsl_rng *rng = NULL;
    if(rng == NULL){
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_default);
    }
    for(p=Np3sub; p<Np3bb; ++p){
        double xrand, yrand, zrand;
        do{
            xrand = xbbmin + gsl_rng_uniform(rng) * (xbbmax-xbbmin);
            yrand = ybbmin + gsl_rng_uniform(rng) * (ybbmax-ybbmin);
            zrand = zbbmin + gsl_rng_uniform(rng) * (zbbmax-zbbmin);
        }while(xrand>=xmin && xrand<xmax && yrand>=ymin && yrand<ymax && zrand>=zmin && zrand<zmax);
        (*xbb)[p] = xrand;
        (*ybb)[p] = yrand;
        (*zbb)[p] = zrand;
    }

    fprintf(stderr, "subbox() %.3fs to pick out %lld(%lld)/%lld=%.2f%%(%.2f%%) particles\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3bb, Np3sub, Np3, 100.*Np3bb/Np3, 100.*Np3sub/Np3);
    return Np3bb;
}

long long int subsphere(double *x, double *y, double *z, long long int Np3, double XYZRL[5], double bbox[6],
        double **xbb, double **ybb, double **zbb, long long int Np3sub_max)
{
    double X = XYZRL[0], Y = XYZRL[1], Z = XYZRL[2], R = XYZRL[3], L=XYZRL[4];
    double R2 = R*R, Vsub = 4*M_PI/3 * R2*R;
    assert(2*R < L);

    double xbbmin = bbox[0], xbbmax = bbox[1];
    double ybbmin = bbox[2], ybbmax = bbox[3];
    double zbbmin = bbox[4], zbbmax = bbox[5];
    double Vbb = (xbbmax-xbbmin) * (ybbmax-ybbmin) * (zbbmax-zbbmin);
    assert(X>=xbbmin && X<=xbbmax);
    assert(Y>=ybbmin && Y<=ybbmax);
    assert(Z>=zbbmin && Z<=zbbmax);

    *xbb = (double *)malloc(Np3sub_max * sizeof(double)); assert(*xbb!=NULL);
    *ybb = (double *)malloc(Np3sub_max * sizeof(double)); assert(*ybb!=NULL);
    *zbb = (double *)malloc(Np3sub_max * sizeof(double)); assert(*zbb!=NULL);
    long long int p, Np3sub;
    clock_t t = clock();
    for(p=0, Np3sub=0; p<Np3 && Np3sub<Np3sub_max; ++p)
        if(pdist2(x[p], y[p], z[p], X, Y, Z, L) < R2){
            (*xbb)[Np3sub] = x[p];
            (*ybb)[Np3sub] = y[p];
            (*zbb)[Np3sub] = z[p];
            ++ Np3sub;
        }
    assert(p == Np3);

    long long int Np3bb = round(Vbb/Vsub * Np3sub);
    *xbb = (double *)realloc(*xbb, Np3bb * sizeof(double)); assert(*xbb!=NULL);
    *ybb = (double *)realloc(*ybb, Np3bb * sizeof(double)); assert(*ybb!=NULL);
    *zbb = (double *)realloc(*zbb, Np3bb * sizeof(double)); assert(*zbb!=NULL);
    static gsl_rng *rng = NULL;
    if(rng == NULL){
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_default);
    }
    for(p=Np3sub; p<Np3bb; ++p){
        double xrand, yrand, zrand;
        do{
            xrand = xbbmin + gsl_rng_uniform(rng) * (xbbmax-xbbmin);
            yrand = ybbmin + gsl_rng_uniform(rng) * (ybbmax-ybbmin);
            zrand = zbbmin + gsl_rng_uniform(rng) * (zbbmax-zbbmin);
        }while(pdist2(xrand, yrand, zrand, X, Y, Z, L) < R2);
        (*xbb)[p] = xrand;
        (*ybb)[p] = yrand;
        (*zbb)[p] = zrand;
    }

    fprintf(stderr, "subsphere() %.3fs to pick out %lld(%lld)/%lld=%.2f%%(%.2f%%) particles\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3bb, Np3sub, Np3, 100.*Np3bb/Np3, 100.*Np3sub/Np3);
    return Np3bb;
}
