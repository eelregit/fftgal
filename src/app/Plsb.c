/* measure P_l(k) from sub-boxes */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "../fft.h"
#include "../gal.h"
#include "../power.h"
#include "../io/io.h"


double H(double a) {
    double Om = 0.29;
    double OL = 1 - Om;
    return 100 * sqrt(Om/(a*a*a) + OL);
}


int main(int argc, char *argv[]) {
    if (argc != 11) {
        fprintf(stderr, "Usage: %s Ng L Nsub wisdom alpha dK indir a id outdir\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    int Ng = atoi(argv[1]); assert(Ng>1 && Ng<=1024);
    double L = atof(argv[2]); assert(L>0 && L<1e4);
    int Nsub = atoi(argv[3]); assert(Nsub>=2 && Nsub<=8);
    char *wisdom = argv[4];
    double alpha = atof(argv[5]); assert(alpha>=0.01 && alpha<=0.1);
    double dK = atof(argv[6]); assert(dK>0 && dK <.1);
    char *indir = argv[7];
    double a = atof(argv[8]); assert(a>0 && a<1.1);
    int id = atoi(argv[9]); assert(id>=1 && id<=1000);
    char *outdir = argv[10];
    int Ngsub = Ng / Nsub; assert(Ng%Nsub == 0);

    fft_t *grid = fft_init(Ng, L, wisdom);

    const int maxlen = 1024;
    char infile[maxlen], outfile[maxlen];
    int retval = snprintf(infile, maxlen, "%s/a%.4f_%04d.mock", indir, a, id);
    assert(retval>=0 && retval<maxlen);
    gal_t *part = gal_loadqpm_cubic(infile, L);

    fprintf(stderr, "\n################## has rsd ##################\n\n");
    gal_t *partrsd = NULL, *partsub = NULL;
    for (int ilos = 0; ilos < 2; ++ilos)  /* skip {0,0,0}, {1,1,1} */
    for (int jlos = 0; jlos < 2; ++jlos)
    for (int klos = (ilos==0 && jlos==0); klos < 2 - (ilos==1 && jlos==1); ++klos) {
        fprintf(stderr, "================== los div ==================\n");
        double los[3] = {ilos, jlos, klos};
        double aH = a * H(a);
        gal_free(partrsd);
        partrsd = gal_rsd(part, los, aH);

        long Nptot = 0;
        for (int isub = 0; isub < Nsub; ++isub)
        for (int jsub = 0; jsub < Nsub; ++jsub)
        for (int ksub = 0; ksub < Nsub; ++ksub) {
            fprintf(stderr, "------------------ sub div ------------------\n");
            double Lsub = L / Nsub;
            double box[7] = {isub*Lsub, (isub+1)*Lsub,
                             jsub*Lsub, (jsub+1)*Lsub,
                             ksub*Lsub, (ksub+1)*Lsub, L};
            gal_free(partsub);
            partsub = gal_subbox(partrsd, box, alpha);
            Nptot += partsub->Np;

            fft_p2k(grid, partsub);

            P(grid, partsub);

            fft_k2x(grid, 0);
            clock_t t = clock();
            double Ngsub3 = (double)Ngsub*Ngsub*Ngsub;
            for (int i = -Ngsub+1; i < Ngsub; ++i)
            for (int j = -Ngsub+1; j < Ngsub; ++j)
            for (int k = -Ngsub+1; k < Ngsub; ++k) {
                double Winv = Ngsub3 / ((double)(Ngsub-abs(i)) * (Ngsub-abs(j)) * (Ngsub-abs(k)));
                F(grid,(i+Ng)%Ng,(j+Ng)%Ng,(k+Ng)%Ng) *= Winv;
            }
            fprintf(stderr, "main() %.3fs on deconvolution\n", (double)(clock()-t)/CLOCKS_PER_SEC);
            fft_x2k(grid, 0);

            retval = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pl_rsd1_los%d%d%d_sb%d%d%d.txt",
                    outdir, a, id, ilos, jlos, klos, isub, jsub, ksub);
            assert(retval>=0 && retval<maxlen);
            Pl(grid, partsub, los, dK, outfile);
            retval = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pmu_rsd1_los%d%d%d_sb%d%d%d.txt",
                    outdir, a, id, ilos, jlos, klos, isub, jsub, ksub);
            assert(retval>=0 && retval<maxlen);
            Pmu(grid, partsub, los, dK, outfile);
        }
        assert(Nptot == part->Np);
    }

    fprintf(stderr, "\n################ without rsd ################\n\n");
    long Nptot = 0;
    double *fk_copy = NULL;
    for (int isub = 0; isub < Nsub; ++isub)
    for (int jsub = 0; jsub < Nsub; ++jsub)
    for (int ksub = 0; ksub < Nsub; ++ksub) {
        fprintf(stderr, "================== sub div ==================\n");
        double Lsub = L / Nsub;
        double box[7] = {isub*Lsub, (isub+1)*Lsub,
                         jsub*Lsub, (jsub+1)*Lsub,
                         ksub*Lsub, (ksub+1)*Lsub, L};
        gal_free(partsub);
        partsub = gal_subbox(part, box, alpha);
        Nptot += partsub->Np;

        fft_p2k(grid, partsub);
        free(fk_copy);
        fk_copy = fft_exportf(grid);

        for (int ilos = 0; ilos < 2; ++ilos)  /* skip {0,0,0}, {1,1,1} */
        for (int jlos = 0; jlos < 2; ++jlos)
        for (int klos = (ilos==0 && jlos==0); klos < 2 - (ilos==1 && jlos==1); ++klos) {
            fprintf(stderr, "------------------ los div ------------------\n");
            fft_importf(grid, fk_copy);
            P(grid, partsub);

            fft_k2x(grid, 0);
            clock_t t = clock();
            double Ngsub3 = (double)Ngsub*Ngsub*Ngsub;
            for (int i = -Ngsub+1; i < Ngsub; ++i)
            for (int j = -Ngsub+1; j < Ngsub; ++j)
            for (int k = -Ngsub+1; k < Ngsub; ++k) {
                double Winv = Ngsub3 / ((double)(Ngsub-abs(i)) * (Ngsub-abs(j)) * (Ngsub-abs(k)));
                F(grid,(i+Ng)%Ng,(j+Ng)%Ng,(k+Ng)%Ng) *= Winv;
            }
            fprintf(stderr, "main() %.3fs on deconvolution\n", (double)(clock()-t)/CLOCKS_PER_SEC);
            fft_x2k(grid, 0);

            double los[3] = {ilos, jlos, klos};
            retval = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pl_rsd0_los%d%d%d_sb%d%d%d.txt",
                    outdir, a, id, ilos, jlos, klos, isub, jsub, ksub);
            assert(retval>=0 && retval<maxlen);
            Pl(grid, partsub, los, dK, outfile);
            retval = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pmu_rsd0_los%d%d%d_sb%d%d%d.txt",
                    outdir, a, id, ilos, jlos, klos, isub, jsub, ksub);
            assert(retval>=0 && retval<maxlen);
            Pmu(grid, partsub, los, dK, outfile);
        }
    }
    assert(Nptot == part->Np);

    free(fk_copy);
    fft_free(grid);
    gal_free(part); gal_free(partrsd); gal_free(partsub);
    return 0;
}
