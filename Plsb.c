/* measure P_l(k) from subboxes */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "fftgal.h"
#include "power.h"
#include "box.h"
#include "io/qpm.h"


double H(double a)
{
    double Om = 0.31;
    double OL = 1 - Om;
    return 100 * sqrt(Om/(a*a*a) + OL);
}


int main(int argc, char *argv[])
{
    if(argc!=7){
        fprintf(stderr, "Usage: %s Ng L wisdom Nsb catalog Np3\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    int Ng = atoi(argv[1]); assert(Ng>1 && Ng<=1024);
    double L = atof(argv[2]); assert(L>0. && L<1e4);
    char *wisdom = argv[3];
    int Nsb = atoi(argv[4]); assert(Nsb>0 && Nsb<=8);
    char *catalog = argv[5];
    int Np3 = atoi(argv[6]); assert(Np3>0 && Np3<=(1<<24));

    fftgal_t *fg = fftgal_init(Ng, L, wisdom);

    double *x, *y, *z, *vx, *vy, *vz, *M;
    int *issat;
    qpm_cubic_mocks_read(catalog, Np3, x, y, z, vx, vy, vz, M, issat);

    double *xd, *yd, *zd;
    double los[3] = {0., 0., 1.};
    double a;
    int catid;
    int ret = sscanf(catalog,
            "/project/projectdirs/boss/galaxy/QPM/dr12d_cubic_mocks/a%lf_%04d.mock",
            &a, &catid); assert(ret==2);
    double aH = a * H(a);
    rsd(x, y, z, vx, vy, vz, Np3, xd, yd, zd, los, aH);
    pbc(xd, yd, zd, Np3, L);

    double *xsb, *ysb, *zsb;
    double xyzlim[6] = {0., 2560., 0., 2560., 0., 2560.};
    int Np3sb = subbox(xd, yd, zd, Np3, xyzlim, xsb, ysb, zsb, Np3);
    assert(Np3sb==Np3);

    fftgal_x2fk(fg, xsb, ysb, zsb, Np3);

    double dK = 0.02;
    Pl(fg, dK, los, "testout.txt");

    return 0;
}
