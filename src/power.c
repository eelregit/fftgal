#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include "fft.h"
#include "gal.h"
#include "power.h"


int Pl(fft_t *grid, gal_t *part, double dK, double los[3], char *file) {
    int Ng = grid->Ng;
    hat(los);
    double KF = 2*M_PI / grid->L;
    double dKinv = KF / dK;
    int Nb = (int)floor(M_SQRT3*(Ng/2) * dKinv) + 1;
    double K[Nb], P0[Nb], P2[Nb], P4[Nb], P6[Nb];
    long N[Nb];
    for (int b = 0; b < Nb; ++b) {
        K[b] = P0[b] = P2[b] = P4[b] = P6[b] = 0;
        N[b] = 0;
    }
    clock_t t = clock();
    double Kval[Ng];
    Kval[0] = 0;
    for (int i = 1; i <= Ng/2; ++i) {
        Kval[i] = i;
        Kval[Ng-i] = - i;
    }
    for (int i = 0; i < Ng; ++i)  /* skip {0,0,0} */
    for (int j = 0; j < Ng; ++j)
    for (int k = (i==0 && j==0); k <= Ng/2; ++k) {
        double Kamp = sqrt(gsl_pow_2(Kval[i]) + gsl_pow_2(Kval[j]) + gsl_pow_2(Kval[k]));
        int b = (int)floor(Kamp * dKinv);
        int count = 1 + (2*k%Ng > 0);
        double delta2 = gsl_pow_2(F_Re(grid,i,j,k)) + gsl_pow_2(F_Im(grid,i,j,k));
        delta2 *= count;
        double mu2 = gsl_pow_2((Kval[i]*los[0] + Kval[j]*los[1] + Kval[k]*los[2]) / Kamp);
        K[b] += count * Kamp;
        P0[b] += delta2;
        P2[b] += delta2 * (1.5*mu2 - 0.5);
        P4[b] += delta2 * ((4.375*mu2 - 3.75)*mu2 + 0.375);
        P6[b] += delta2 * (((14.4375*mu2 - 19.6875)*mu2 + 6.5625)*mu2 - 0.3125);
        N[b] += count;
    }
    fprintf(stderr, "Pl() %.3fs on binning, los[]={%.3f,%.3f,%.3f}\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, los[0], los[1], los[2]);

    long Ntot = 0;
    double alpha;
    if (part->rand == NULL)
        alpha = 0;
    else
        alpha = (double)part->Np / part->rand->Np;
    for (int b = 0; b < Nb; ++b) {
        K[b] *= KF / N[b];
        P0[b] *= 1 / part->V / N[b];
        P0[b] -= (1 + alpha) * part->V / part->Np;
        P2[b] *= 5 / part->V / N[b];
        P4[b] *= 9 / part->V / N[b];
        P6[b] *= 13 / part->V / N[b];
        Ntot += N[b];
    }
    assert(Ntot == Ng*Ng*Ng-1);

    FILE *fp = fopen(file, "w"); assert(fp!=NULL);
    fprintf(fp, "# Nb %d\n", Nb);
    fprintf(fp, "# dK %g\n", dK);
    fprintf(fp, "# Ng %d\n", grid->Ng);
    fprintf(fp, "# L %g\n", grid->L);
    fprintf(fp, "# Np %ld\n", part->Np);
    fprintf(fp, "# V %g\n", part->V);
    fprintf(fp, "# alpha %g\n", alpha);
    fprintf(fp, "# los[3] %g %g %g\n", los[0], los[1], los[2]);
    fprintf(fp, "# K P0 P2 P4 P6 N\n");
    for (int b = 0; b < Nb; ++b)
        fprintf(fp, "%e %e % e % e % e %ld\n", K[b], P0[b], P2[b], P4[b], P6[b], N[b]);
    fclose(fp);
    fprintf(stderr, "Pl() saved to %s\n", file);

    return Nb;
}
