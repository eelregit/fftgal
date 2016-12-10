#include <math.h>


void rsd(double *x, double *y, double *z, double *vx, double *vy, double *vz,
        long long int Np3, double *xd, double *yd, double *zd,
        double los[3], double aH)
{
    xd = (double *)malloc(Np3 * sizeof(double)); assert(xd!=NULL);
    yd = (double *)malloc(Np3 * sizeof(double)); assert(yd!=NULL);
    zd = (double *)malloc(Np3 * sizeof(double)); assert(zd!=NULL);

    double los_ampl = sqrt(los[1]*los[1], los[2]*los[2], los[3]*los[3]);
    los[1] /= los_ampl;
    los[2] /= los_ampl;
    los[3] /= los_ampl;

    double aHinv = 1 / aH;
    long long int p;
    for(p=0; p<Np3; ++p){
        vdotlos = los[1]*vx[p] + los[2]*vy[p] + los[3]*vz[p];
        xd[p] = x[p] + aHinv * vdotlos * los[1];
        yd[p] = y[p] + aHinv * vdotlos * los[2];
        zd[p] = z[p] + aHinv * vdotlos * los[3];
    }
}
