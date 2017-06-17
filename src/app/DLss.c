/* measure super-sample modes from sub-spheres */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include "../fft.h"
#include "../gal.h"
#include "../io/io.h"


/* use bias!=1 to rescale to matter overdensity if it's known,
 * otherwise output is in tracer overdensity */
const double bias = 1;


static double tophat(double KR) {
    return 3 * (sin(KR) - KR*cos(KR)) / gsl_pow_3(KR);  /* KR != 0 */
}


int main(int argc, char *argv[]) {
    if (argc != 9) {
        fprintf(stderr, "Usage: %s Ng L Nsub wisdom indir a id outdir\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    int Ng = atoi(argv[1]); assert(Ng>1 && Ng<=1024);
    double L = atof(argv[2]); assert(L>0 && L<1e4);
    int Nsub = atoi(argv[3]); assert(Nsub>=2 && Nsub<=8);
    char *wisdom = argv[4];
    char *indir = argv[5];
    double a = atof(argv[6]); assert(a>0 && a<1.1);
    int id = atoi(argv[7]); assert(id>=1 && id<=1000);
    char *outdir = argv[8];
    int Ngsub = Ng / Nsub; assert(Ng%Nsub == 0);

    fft_t *grid = fft_init(Ng, L, wisdom);

    const int maxlen = 1024;
    char infile[maxlen];
    int retval = snprintf(infile, maxlen, "%s/a%.4f_%04d.mock", indir, a, id);
    assert(retval>=0 && retval<maxlen);
    gal_t *part = gal_loadqpm_cubic(infile, L);

    double *fk_copy = NULL;
    for (int ilos = 0; ilos < 2; ++ilos)  /* skip {0,0,0}, {1,1,1} */
    for (int jlos = 0; jlos < 2; ++jlos)
    for (int klos = (ilos==0 && jlos==0); klos < 2 - (ilos==1 && jlos==1); ++klos) {
        fprintf(stderr, "================== los div ==================\n");
        double los[3] = {ilos, jlos, klos};
        hat(los);

        double DeltaL[3][Nsub*Nsub*Nsub];
        for (int iL = 0; iL <= 2; ++iL) {
            if (fk_copy == NULL) {
                double offset[3] = {0, 0, 0};
                fft_p2g(grid, part, offset);
                fft_x2k(grid, 0);
                fk_copy = fft_exportf(grid);
            } else {
                fft_importf(grid, fk_copy);
            }

            double Kval[Ng];
            Kval[0] = 0;
            for (int i = 1; i <= Ng/2; ++i) {
                Kval[i] = i;
                Kval[Ng-i] = - i;
            }
            grid->f[0] = 0;
            for (int i = 0; i < Ng; ++i)  /* skip {0,0,0} */
            for (int j = 0; j < Ng; ++j)
            for (int k = (i==0 && j==0); k <= Ng/2; ++k) {
                double Kamp = sqrt(gsl_pow_2(Kval[i]) + gsl_pow_2(Kval[j]) + gsl_pow_2(Kval[k]));
                double mu2 = gsl_pow_2((Kval[i]*los[0] + Kval[j]*los[1] + Kval[k]*los[2]) / Kamp);
                double Legendre[3] = {1, 1.5*mu2 - 0.5, (4.375*mu2 - 3.75)*mu2 + 0.375};
                double W = tophat(Kamp * M_PI*M_SQRT3 / Nsub);
                F_Re(grid,i,j,k) *= Legendre[iL] * W;
                F_Im(grid,i,j,k) *= Legendre[iL] * W;
            }

            fft_k2x(grid, 0);
            for (int isub = 0; isub < Nsub; ++isub)
            for (int jsub = 0; jsub < Nsub; ++jsub)
            for (int ksub = 0; ksub < Nsub; ++ksub)
                DeltaL[iL][(isub*Nsub + jsub)*Nsub + ksub]
                    = F(grid, isub*Ngsub, jsub*Ngsub, ksub*Ngsub) / bias;
        }

        char outfile[maxlen];
        retval = snprintf(outfile, maxlen, "%s/a%.4f_%04d/DLss_los%d%d%d.txt",
                outdir, a, id, ilos, jlos, klos);
        assert(retval>=0 && retval<maxlen);
        FILE *fp = fopen(outfile, "w"); assert(fp!=NULL);
        fprintf(fp, "# ijk_sub Delta_0 Delta_2 Delta_4\n");
        for (int isub = 0; isub < Nsub; ++isub)
        for (int jsub = 0; jsub < Nsub; ++jsub)
        for (int ksub = 0; ksub < Nsub; ++ksub)
            fprintf(fp, "%d%d%d % e % e % e\n", isub, jsub, ksub,
                    DeltaL[0][(isub*Nsub + jsub)*Nsub + ksub],
                    DeltaL[1][(isub*Nsub + jsub)*Nsub + ksub],
                    DeltaL[2][(isub*Nsub + jsub)*Nsub + ksub]);
        fclose(fp);
    }

    free(fk_copy);
    fft_free(grid);
    gal_free(part);
    return 0;
}
