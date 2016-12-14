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
    if(argc!=10){
        fprintf(stderr, "Usage: %s Ng L wisdom Nsb catdir a catid Np3 outdir\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    int Ng = atoi(argv[1]); assert(Ng>1 && Ng<=1024);
    double L = atof(argv[2]); assert(L>0. && L<1e4);
    char *wisdom = argv[3];
    int Nsb = atoi(argv[4]); assert(Nsb>0 && Nsb<=8);
    char *catdir = argv[5];
    double a = atof(argv[6]); assert(a>0. && a<1.1);
    int catid = atoi(argv[7]); assert(catid>=1 && catid<=1000);
    int Np3 = atoi(argv[8]); assert(Np3>0 && Np3<=(1<<24));
    char *outdir = argv[9];

    const int maxlen = 1024;
    char catalog[maxlen];
    int ret = snprintf(catalog, maxlen, "%s/a%.4f_%04d.mock", catdir, a, catid);
    assert(ret>=0 && ret<maxlen);
    double *x, *y, *z, *vx, *vy, *vz, *M;
    int *issat;
    qpm_cubic_mocks_read(catalog, Np3, &x, &y, &z, &vx, &vy, &vz, &M, &issat);

    double Lsb = L / Nsb;
    fftgal_t *fg = fftgal_init(Ng, Lsb, wisdom);

    double *xd=NULL, *yd=NULL, *zd=NULL;
    double *xsb=NULL, *ysb=NULL, *zsb=NULL;
    for(int ilos=0; ilos<2; ++ilos)
    for(int jlos=0; jlos<2; ++jlos)
    for(int klos=0; klos<2-(ilos==1 && jlos==1); ++klos){ /* skip los[] = {1,1,1} */
        fprintf(stderr, "================== los div ==================\n");
        if(xd!=NULL) free(xd);
        if(yd!=NULL) free(yd);
        if(zd!=NULL) free(zd);
        double los[3] = {ilos, jlos, klos};
        double aH = a * H(a);
        rsd(x, y, z, vx, vy, vz, Np3, &xd, &yd, &zd, los, aH);
        pbc(xd, yd, zd, Np3, L);

        for(int isb=0; isb<Nsb; ++isb)
        for(int jsb=0; jsb<Nsb; ++jsb)
        for(int ksb=0; ksb<Nsb; ++ksb){
            fprintf(stderr, "---------------- sub-box div ----------------\n");
            if(xsb!=NULL) free(xsb);
            if(ysb!=NULL) free(ysb);
            if(zsb!=NULL) free(zsb);
            double xyzlim[6] = {isb*Lsb, (isb+1)*Lsb,
                                jsb*Lsb, (jsb+1)*Lsb,
                                ksb*Lsb, (ksb+1)*Lsb};
            int Np3sb = subbox(xd, yd, zd, Np3, xyzlim, &xsb, &ysb, &zsb, Np3);

            fftgal_x2fk(fg, xsb, ysb, zsb, Np3sb);

            char outfile[maxlen];
            ret = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pl_los%d%d%d_sb%d%d%d.txt",
                    outdir, a, catid, ilos, jlos, klos, isb, jsb, ksb);
            assert(ret>=0 && ret<maxlen);
            double dK = 0.02;
            Pl(fg, dK, los, outfile);
        }
    }

    free(x); free(y); free(z); free(vx); free(vy); free(vz); free(M); free(issat);
    free(xd); free(yd); free(zd);
    free(xsb); free(ysb); free(zsb);
    fftgal_kill(fg);
    return 0;
}
