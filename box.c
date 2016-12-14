#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "box.h"


void pbc(double *x, double *y, double *z, long long int Np3, double L)
{
    double Linv = 1 / L;
    clock_t t = clock();
    for(long long int p=0; p<Np3; ++p){
        x[p] -= L * floor(x[p] * Linv);
        y[p] -= L * floor(y[p] * Linv);
        z[p] -= L * floor(z[p] * Linv);
    }
    fprintf(stderr, "pbc() %.3f sec to wrap around\n", (double)(clock()-t)/CLOCKS_PER_SEC);
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
    fprintf(stderr, "subbox() %.3f to pick out %lld/%lld particles (%.2f%%)\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3sb, Np3, percentage);
    if(Np3sb!=Np3 && percentage>99.)
        fprintf(stderr, "warning: subbox() thinks something is leaking\n");
    return Np3sb;
}
