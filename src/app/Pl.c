/* measure P_l(k) from one periodic box */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../fftgal.h"
#include "../power.h"
#include "../geom.h"
#include "../io/qpm.h"


double H(double a)
{
    double Om = 0.29;
    double OL = 1 - Om;
    return 100 * sqrt(Om/(a*a*a) + OL);
}


int main(int argc, char *argv[])
{
    if(argc!=10){
        fprintf(stderr, "Usage: %s Ng L dK wisdom fold catdir a catid outdir\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    int Ng = atoi(argv[1]); assert(Ng>1 && Ng<=1024);
    double L = atof(argv[2]); assert(L>0. && L<1e4);
    double dK = atof(argv[3]); assert(dK>0. && dK <.1);
    char *wisdom = argv[4];
    int fold = atoi(argv[5]); assert(fold>0 && fold<=8);
    char *catdir = argv[6];
    double a = atof(argv[7]); assert(a>0. && a<1.1);
    int catid = atoi(argv[8]); assert(catid>=1 && catid<=1000);
    char *outdir = argv[9];

    const int maxlen = 1024;
    char catalog[maxlen];
    int ret = snprintf(catalog, maxlen, "%s/a%.4f_%04d.mock", catdir, a, catid);
    assert(ret>=0 && ret<maxlen);
    double *x, *y, *z, *vx, *vy, *vz, *M;
    int *issat;
    int Np3 = qpm_cubic_mocks_load(catalog, &x, &y, &z, &vx, &vy, &vz, &M, &issat);
    assert(Np3>0);
    pbc(x, y, z, Np3, L);

    fftgal_t *fg = fftgal_init(Ng, L, -1, fold, wisdom);

    fprintf(stderr, "\n################## has rsd ##################\n\n");
    double *xd=NULL, *yd=NULL, *zd=NULL;
    for(int ilos=0; ilos<2; ++ilos)
    for(int jlos=0; jlos<2; ++jlos)
    for(int klos=0+(ilos==0 && jlos==0); klos<2-(ilos==1 && jlos==1); ++klos){ /* skip {0,0,0}, {1,1,1} */
        fprintf(stderr, "================== los div ==================\n");
        if(xd!=NULL) free(xd);
        if(yd!=NULL) free(yd);
        if(zd!=NULL) free(zd);
        double los[3] = {ilos, jlos, klos};
        double aH = a * H(a);
        rsd(x, y, z, vx, vy, vz, Np3, &xd, &yd, &zd, los, aH);
        pbc(xd, yd, zd, Np3, L);

        fftgal_x2fk(fg, xd, yd, zd, Np3);

        char outfile[maxlen];
        ret = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pl_rsd1_los%d%d%d.txt",
                outdir, a, catid, ilos, jlos, klos);
        assert(ret>=0 && ret<maxlen);
        Pl(fg, dK, los, outfile);
    }

    fprintf(stderr, "\n################ without rsd ################\n\n");
    fftgal_x2fk(fg, x, y, z, Np3);

    for(int ilos=0; ilos<2; ++ilos)
    for(int jlos=0; jlos<2; ++jlos)
    for(int klos=0+(ilos==0 && jlos==0); klos<2-(ilos==1 && jlos==1); ++klos){ /* skip {0,0,0}, {1,1,1} */
        fprintf(stderr, "------------------ los div ------------------\n");
        double los[3] = {ilos, jlos, klos};

        char outfile[maxlen];
        ret = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pl_rsd0_los%d%d%d.txt",
                outdir, a, catid, ilos, jlos, klos);
        assert(ret>=0 && ret<maxlen);
        Pl(fg, dK, los, outfile);
    }

    free(x); free(y); free(z); free(vx); free(vy); free(vz); free(M); free(issat);
    free(xd); free(yd); free(zd);
    fftgal_free(fg);
    return 0;
}
