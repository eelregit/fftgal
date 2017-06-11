/* measure P_l(k) from one periodic box */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
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
    if (argc != 10) {
        fprintf(stderr, "Usage: %s Ng L fold wisdom dK indir a id outdir\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    int Ng = atoi(argv[1]); assert(Ng>1 && Ng<=1024);
    double L = atof(argv[2]); assert(L>0 && L<1e4);
    int fold = atoi(argv[3]); assert(fold>0 && fold<=8);
    char *wisdom = argv[4];
    double dK = atof(argv[5]); assert(dK>0 && dK <.1);
    char *indir = argv[6];
    double a = atof(argv[7]); assert(a>0 && a<1.1);
    int id = atoi(argv[8]); assert(id>=1 && id<=1000);
    char *outdir = argv[9];

    double Lfd = L / fold;
    fft_t *grid = fft_init(Ng, Lfd, wisdom);

    const int maxlen = 1024;
    char infile[maxlen];
    int retval = snprintf(infile, maxlen, "%s/a%.4f_%04d.mock", indir, a, id);
    assert(retval>=0 && retval<maxlen);
    gal_t *part = gal_loadqpm_cubic(infile, L);

    fprintf(stderr, "\n################## has rsd ##################\n\n");
    gal_t *partrsd = NULL;
    for (int ilos = 0; ilos < 2; ++ilos)  /* skip {0,0,0}, {1,1,1} */
    for (int jlos = 0; jlos < 2; ++jlos)
    for (int klos = (ilos==0 && jlos==0); klos < 2 - (ilos==1 && jlos==1); ++klos) {
        fprintf(stderr, "================== los div ==================\n");
        double los[3] = {ilos, jlos, klos};
        double aH = a * H(a);
        gal_free(partrsd);
        partrsd = gal_rsd(part, los, aH);

        fft_p2k(grid, partrsd);

        char outfile[maxlen];
        retval = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pl_rsd1_los%d%d%d_fd%d.txt",
                outdir, a, id, ilos, jlos, klos, fold);
        assert(retval>=0 && retval<maxlen);
        Pl(grid, partrsd, los, dK, outfile);
    }

    fprintf(stderr, "\n################ without rsd ################\n\n");

    fft_p2k(grid, part);

    for (int ilos = 0; ilos < 2; ++ilos)  /* skip {0,0,0}, {1,1,1} */
    for (int jlos = 0; jlos < 2; ++jlos)
    for (int klos = (ilos==0 && jlos==0); klos < 2 - (ilos==1 && jlos==1); ++klos) {
        fprintf(stderr, "================== los div ==================\n");
        double los[3] = {ilos, jlos, klos};
        char outfile[maxlen];
        retval = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Pl_rsd0_los%d%d%d_fd%d.txt",
                outdir, a, id, ilos, jlos, klos, fold);
        assert(retval>=0 && retval<maxlen);
        Pl(grid, part, los, dK, outfile);
    }

    fft_free(grid);
    gal_free(part); gal_free(partrsd);
    return 0;
}
