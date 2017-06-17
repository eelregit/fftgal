#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include "fft.h"
#include "gal.h"
#include "power.h"


void P(fft_t *grid, gal_t *part) {
    clock_t t = clock();
    int Ng = grid->Ng;
    double Vinv = 1 / part->V;
    double alpha = (part->rand == NULL) ? 0 : (double)part->Np / part->rand->Np;
    double Pshot = (1 + alpha) * part->V / part->Np;  /* subtract shot noise */
    F_Re(grid,0,0,0) = F_Im(grid,0,0,0) = 0;
    for (int i = 0; i < Ng; ++i)  /* skip {0,0,0} */
    for (int j = 0; j < Ng; ++j)
    for (int k = (i==0 && j==0); k <= Ng/2; ++k) {
        double delta2 = gsl_pow_2(F_Re(grid,i,j,k)) + gsl_pow_2(F_Im(grid,i,j,k));
        F_Re(grid,i,j,k) = delta2 * Vinv - Pshot;
        F_Im(grid,i,j,k) = 0;
    }
    fprintf(stderr, "P() %.3fs\n", (double)(clock()-t)/CLOCKS_PER_SEC);
}


void Pl(fft_t *grid, gal_t *part, double los[3], double dK, char *file) {
    clock_t t = clock();
    int Ng = grid->Ng;
    hat(los);
    double KF = 2*M_PI / grid->L;
    double dKinv = KF / dK;
    int NK = floor(M_SQRT3*(Ng/2) * dKinv) + 1;
    int Nl = 4;
    double K[NK], P[NK][Nl];
    long N[NK];
    for (int iK = 0; iK < NK; ++iK) {
        K[iK] = 0;
        N[iK] = 0;
        for (int il = 0; il < Nl; ++il)
            P[iK][il] = 0;
    }
    double Kval[Ng];
    Kval[0] = 0;
    for (int i = 1; i <= Ng/2; ++i) {
        Kval[i] = i;
        Kval[Ng-i] = - i;
    }
    for (int i = 0; i < Ng; ++i)  /* skip {0,0,0} */
    for (int j = 0; j < Ng; ++j)
    for (int k = (i==0 && j==0); k <= Ng/2; ++k) {
        int count = 1 + (2*k%Ng > 0);
        double delta2 = F_Re(grid,i,j,k);
        delta2 *= count;
        double Kamp = sqrt(gsl_pow_2(Kval[i]) + gsl_pow_2(Kval[j]) + gsl_pow_2(Kval[k]));
        int iK = floor(Kamp * dKinv);
        double mu2 = gsl_pow_2((Kval[i]*los[0] + Kval[j]*los[1] + Kval[k]*los[2]) / Kamp);
        K[iK] += count * Kamp;
        N[iK] += count;
        P[iK][0] += delta2;
        P[iK][1] += delta2 * (1.5*mu2 - 0.5);
        P[iK][2] += delta2 * ((4.375*mu2 - 3.75)*mu2 + 0.375);
        P[iK][3] += delta2 * (((14.4375*mu2 - 19.6875)*mu2 + 6.5625)*mu2 - 0.3125);
    }

    long Ntot = 0;
    for (int iK = 0; iK < NK; ++iK) {
        K[iK] *= KF / N[iK];
        for (int il = 0; il < Nl; ++il)
            P[iK][il] *= (4*il + 1.0) / N[iK];
        Ntot += N[iK];
    }
    assert(Ntot == Ng*Ng*Ng-1);

    double alpha = (part->rand == NULL) ? 0 : (double)part->Np / part->rand->Np;
    FILE *fp = fopen(file, "w"); assert(fp!=NULL);
    fprintf(fp, "# NK %d\n", NK);
    fprintf(fp, "# dK %g\n", dK);
    fprintf(fp, "# Nl %d\n", Nl);
    fprintf(fp, "# Ng %d\n", grid->Ng);
    fprintf(fp, "# L %g\n", grid->L);
    fprintf(fp, "# Np %ld\n", part->Np);
    fprintf(fp, "# V %g\n", part->V);
    fprintf(fp, "# alpha %g\n", alpha);
    fprintf(fp, "# los[3] %g %g %g\n", los[0], los[1], los[2]);
    fprintf(fp, "# K N");
    for (int il = 0; il < Nl; ++il)
        fprintf(fp, " P%d", 2*il);
    fprintf(fp, "\n");
    for (int iK = 0; iK < NK; ++iK) {
        fprintf(fp, "%e %9ld", K[iK], N[iK]);
        for (int il = 0; il < Nl; ++il)
            fprintf(fp, " % e", P[iK][il]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    fprintf(stderr, "Pl() %.3fs, saved to %s\n", (double)(clock()-t)/CLOCKS_PER_SEC, file);
}


void Pmu(fft_t *grid, gal_t *part, double los[3], double dK, char *file) {
    clock_t t = clock();
    int Ng = grid->Ng;
    hat(los);
    double KF = 2*M_PI / grid->L;
    double dKinv = KF / dK;
    int NK = floor(M_SQRT3*(Ng/2) * dKinv) + 1;
    int Nmu = 3;
    double K[NK][Nmu+1], P[NK][Nmu+1];
    long N[NK][Nmu+1];
    for (int iK = 0; iK < NK; ++iK)
        for (int imu = 0; imu <= Nmu; ++imu) {
            K[iK][imu] = P[iK][imu] = 0;
            N[iK][imu] = 0;
        }
    double Kval[Ng];
    Kval[0] = 0;
    for (int i = 1; i <= Ng/2; ++i) {
        Kval[i] = i;
        Kval[Ng-i] = - i;
    }
    for (int i = 0; i < Ng; ++i)  /* skip {0,0,0} */
    for (int j = 0; j < Ng; ++j)
    for (int k = (i==0 && j==0); k <= Ng/2; ++k) {
        int count = 1 + (2*k%Ng > 0);
        double delta2 = F_Re(grid,i,j,k);
        delta2 *= count;
        double Kamp = sqrt(gsl_pow_2(Kval[i]) + gsl_pow_2(Kval[j]) + gsl_pow_2(Kval[k]));
        int iK = floor(Kamp * dKinv);
        double mu = (Kval[i]*los[0] + Kval[j]*los[1] + Kval[k]*los[2]) / Kamp;
        int imu = floor(fabs(mu * Nmu));
        K[iK][imu] += count * Kamp;
        N[iK][imu] += count;
        P[iK][imu] += delta2;
    }
    for (int iK = 0; iK < NK; ++iK) {
        K[iK][Nmu-1] += K[iK][Nmu];  /* corner case mu = 1 */
        N[iK][Nmu-1] += N[iK][Nmu];
        P[iK][Nmu-1] += P[iK][Nmu];
    }

    long Ntot = 0;
    for (int iK = 0; iK < NK; ++iK)
        for (int imu = 0; imu < Nmu; ++imu)
            if (N[iK][imu]) {
                K[iK][imu] *= KF / N[iK][imu];
                P[iK][imu] /= N[iK][imu];
                Ntot += N[iK][imu];
            }
    assert(Ntot == Ng*Ng*Ng-1);

    double alpha = (part->rand == NULL) ? 0 : (double)part->Np / part->rand->Np;
    FILE *fp = fopen(file, "w"); assert(fp!=NULL);
    fprintf(fp, "# NK %d\n", NK);
    fprintf(fp, "# dK %g\n", dK);
    fprintf(fp, "# Nmu %d\n", Nmu);
    fprintf(fp, "# Ng %d\n", grid->Ng);
    fprintf(fp, "# L %g\n", grid->L);
    fprintf(fp, "# Np %ld\n", part->Np);
    fprintf(fp, "# V %g\n", part->V);
    fprintf(fp, "# alpha %g\n", alpha);
    fprintf(fp, "# los[3] %g %g %g\n", los[0], los[1], los[2]);
    fprintf(fp, "#");
    for (int imu = 0; imu < Nmu; ++imu)
        fprintf(fp, " K%d N%d P%d", imu, imu, imu);
    fprintf(fp, "\n");
    for (int iK = 0; iK < NK; ++iK) {
        fprintf(fp, "%e %9ld % e", K[iK][0], N[iK][0], P[iK][0]);
        for (int imu = 1; imu < Nmu; ++imu)
            fprintf(fp, " %e %9ld % e", K[iK][imu], N[iK][imu], P[iK][imu]);
        fprintf(fp, "\n");
    }
    fclose(fp);
    fprintf(stderr, "Pmu() %.3fs, saved to %s\n", (double)(clock()-t)/CLOCKS_PER_SEC, file);
}
