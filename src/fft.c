#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include "gal.h"
#include "fft.h"


fft_t *fft_init(int Ng, double L, char *wisdom) {
    fft_t *self = (fft_t *)malloc(sizeof(fft_t)); assert(self!=NULL);
    for (int i = 0; i < 3; ++i)
        self->offset[i] = 0;
    self->Ng = Ng;
    self->L = L;
    fprintf(stderr, "fft_init() making a %d^3 grid in a %g^3 box\n", Ng, L);

    int Ng_pad = 2 * (Ng / 2 + 1);
    self->Ng_pad = Ng_pad;
    long Ng3_pad = (long)Ng * Ng * Ng_pad;
    self->Ng3_pad  = Ng3_pad;
    double *fx = (double *)fftw_malloc(Ng3_pad * sizeof(double)); assert(fx!=NULL);
    fftw_complex *fk = (fftw_complex *)fx;
    self->f = fx;

    if (!strcmp(wisdom, "FFTW_ESTIMATE")) {
        fprintf(stderr, "fft_init() using FFTW_ESTIMATE\n");
        self->x2k = fftw_plan_dft_r2c_3d(Ng, Ng, Ng, fx, fk, FFTW_ESTIMATE);
        self->k2x = fftw_plan_dft_c2r_3d(Ng, Ng, Ng, fk, fx, FFTW_ESTIMATE);
    } else {
        int retval = fftw_import_wisdom_from_filename(wisdom);
        if (retval) {
            fprintf(stderr, "fft_init() imported wisdom from %s\n", wisdom);
            self->x2k = fftw_plan_dft_r2c_3d(Ng, Ng, Ng, fx, fk, FFTW_ESTIMATE);
            self->k2x = fftw_plan_dft_c2r_3d(Ng, Ng, Ng, fk, fx, FFTW_ESTIMATE);
        } else {
            clock_t t = clock();
            self->x2k = fftw_plan_dft_r2c_3d(Ng, Ng, Ng, fx, fk, FFTW_MEASURE);
            self->k2x = fftw_plan_dft_c2r_3d(Ng, Ng, Ng, fk, fx, FFTW_MEASURE);
            retval = fftw_export_wisdom_to_filename(wisdom); assert(retval);
            fprintf(stderr, "fft_init() exported wisdom to %s\n", wisdom);
            fprintf(stderr, "fft_init() %.3fs on FFTW_MEASURE\n",
                    (double)(clock()-t)/CLOCKS_PER_SEC);
        }
    }

    return self;
}


void fft_free(fft_t *self) {
    if (self != NULL)
        fftw_free(self->f);
    free(self);
}


void fft_p2g(fft_t *self, gal_t *part, double offset[3]) {
    int Ng = self->Ng;
    for (int i = 0; i < 3; ++i) {
        offset[i] = remainder(offset[i], Ng);
        self->offset[i] = offset[i];
    }

    clock_t t = clock();
    for (long g = 0; g < self->Ng3_pad; ++g)  /* must zero after fftw_plan */
        self->f[g] = 0;
    double Hinv = Ng / self->L;
    for (long p = 0; p < part->Np; ++p) {
        double x = part->x[p] * Hinv - offset[0];
        double y = part->y[p] * Hinv - offset[1];
        double z = part->z[p] * Hinv - offset[2];
        int gi1 = (int)floor(x);
        int gj1 = (int)floor(y);
        int gk1 = (int)floor(z);
        double dx = x - gi1;
        double dy = y - gj1;
        double dz = z - gk1;
        gi1 = (gi1 % Ng + Ng) % Ng;
        gj1 = (gj1 % Ng + Ng) % Ng;
        gk1 = (gk1 % Ng + Ng) % Ng;
        int gi[4] = {(gi1 - 1 + Ng) % Ng, gi1, (gi1 + 1) % Ng, (gi1 + 2) % Ng};
        int gj[4] = {(gj1 - 1 + Ng) % Ng, gj1, (gj1 + 1) % Ng, (gj1 + 2) % Ng};
        int gk[4] = {(gk1 - 1 + Ng) % Ng, gk1, (gk1 + 1) % Ng, (gk1 + 2) % Ng};
        double wx[4] = {gsl_pow_3(1 - dx) / 6, ((3*dx - 6)*gsl_pow_2(dx) + 4) / 6,
                    (((-3*dx + 3)*dx + 3)*dx + 1) / 6, gsl_pow_3(dx) / 6};
        double wy[4] = {gsl_pow_3(1 - dy) / 6, ((3*dy - 6)*gsl_pow_2(dy) + 4) / 6,
                    (((-3*dy + 3)*dy + 3)*dy + 1) / 6, gsl_pow_3(dy) / 6};
        double wz[4] = {gsl_pow_3(1 - dz) / 6, ((3*dz - 6)*gsl_pow_2(dz) + 4) / 6,
                    (((-3*dz + 3)*dz + 3)*dz + 1) / 6, gsl_pow_3(dz) / 6};
        for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k)
            F(self, gi[i], gj[j], gk[k]) += wx[i] * wy[j] * wz[k];
    }
    if (part->rand != NULL) {
        gal_t *rand = part->rand;
        double alpha = (double)part->Np / rand->Np;
        for (long p = 0; p < rand->Np; ++p) {
            double x = rand->x[p] * Hinv - offset[0];
            double y = rand->y[p] * Hinv - offset[1];
            double z = rand->z[p] * Hinv - offset[2];
            int gi1 = (int)floor(x);
            int gj1 = (int)floor(y);
            int gk1 = (int)floor(z);
            double dx = x - gi1;
            double dy = y - gj1;
            double dz = z - gk1;
            gi1 = (gi1 % Ng + Ng) % Ng;
            gj1 = (gj1 % Ng + Ng) % Ng;
            gk1 = (gk1 % Ng + Ng) % Ng;
            int gi[4] = {(gi1 - 1 + Ng) % Ng, gi1, (gi1 + 1) % Ng, (gi1 + 2) % Ng};
            int gj[4] = {(gj1 - 1 + Ng) % Ng, gj1, (gj1 + 1) % Ng, (gj1 + 2) % Ng};
            int gk[4] = {(gk1 - 1 + Ng) % Ng, gk1, (gk1 + 1) % Ng, (gk1 + 2) % Ng};
            double wx[4] = {gsl_pow_3(1 - dx) / 6, ((3*dx - 6)*gsl_pow_2(dx) + 4) / 6,
                        (((-3*dx + 3)*dx + 3)*dx + 1) / 6, gsl_pow_3(dx) / 6};
            double wy[4] = {gsl_pow_3(1 - dy) / 6, ((3*dy - 6)*gsl_pow_2(dy) + 4) / 6,
                        (((-3*dy + 3)*dy + 3)*dy + 1) / 6, gsl_pow_3(dy) / 6};
            double wz[4] = {gsl_pow_3(1 - dz) / 6, ((3*dz - 6)*gsl_pow_2(dz) + 4) / 6,
                        (((-3*dz + 3)*dz + 3)*dz + 1) / 6, gsl_pow_3(dz) / 6};
            for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                F(self, gi[i], gj[j], gk[k]) -= alpha * wx[i] * wy[j] * wz[k];
        }
    }
    double H3 = gsl_pow_3(self->L / Ng);
    double nbar = part->Np / part->V;
    double Npcinv = 1 / nbar / H3;  /* inverse of Np/cell */
    for (long g = 0; g < self->Ng3_pad; ++g)
        self->f[g] *= Npcinv;
    fprintf(stderr, "fft_p2g() %.3fs\n", (double)(clock()-t)/CLOCKS_PER_SEC);
}


double *fft_g2p(fft_t *self, gal_t *pos) {
    double *f_interp = (double *)malloc(pos->Np * sizeof(double)); assert(f_interp!=NULL);
    clock_t t = clock();
    int Ng = self->Ng;
    double Hinv = Ng / self->L;
    for (long p = 0; p < pos->Np; ++p) {
        double x = pos->x[p] * Hinv;
        double y = pos->y[p] * Hinv;
        double z = pos->z[p] * Hinv;
        int gi1 = (int)floor(x);
        int gj1 = (int)floor(y);
        int gk1 = (int)floor(z);
        double dx = x - gi1;
        double dy = y - gj1;
        double dz = z - gk1;
        gi1 = (gi1 % Ng + Ng) % Ng;
        gj1 = (gj1 % Ng + Ng) % Ng;
        gk1 = (gk1 % Ng + Ng) % Ng;
        int gi[4] = {(gi1 - 1 + Ng) % Ng, gi1, (gi1 + 1) % Ng, (gi1 + 2) % Ng};
        int gj[4] = {(gj1 - 1 + Ng) % Ng, gj1, (gj1 + 1) % Ng, (gj1 + 2) % Ng};
        int gk[4] = {(gk1 - 1 + Ng) % Ng, gk1, (gk1 + 1) % Ng, (gk1 + 2) % Ng};
        double wx[4] = {gsl_pow_3(1 - dx) / 6, ((3*dx - 6)*gsl_pow_2(dx) + 4) / 6,
                    (((-3*dx + 3)*dx + 3)*dx + 1) / 6, gsl_pow_3(dx) / 6};
        double wy[4] = {gsl_pow_3(1 - dy) / 6, ((3*dy - 6)*gsl_pow_2(dy) + 4) / 6,
                    (((-3*dy + 3)*dy + 3)*dy + 1) / 6, gsl_pow_3(dy) / 6};
        double wz[4] = {gsl_pow_3(1 - dz) / 6, ((3*dz - 6)*gsl_pow_2(dz) + 4) / 6,
                    (((-3*dz + 3)*dz + 3)*dz + 1) / 6, gsl_pow_3(dz) / 6};
        for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k)
            f_interp[p] = F(self, gi[i], gj[j], gk[k]) * wx[i] * wy[j] * wz[k];
    }
    fprintf(stderr, "fft_g2p() %.3fs\n", (double)(clock()-t)/CLOCKS_PER_SEC);
    return f_interp;
}


void fft_x2k(fft_t *self, int offset_phase_on) {
    clock_t t = clock();
    fftw_execute(self->x2k);
    double H3 = gsl_pow_3(self->L / self->Ng);
    for (long g = 0; g < self->Ng3_pad; ++g)
        self->f[g] *= H3;
    if (offset_phase_on)
        fft_offset_phase(self, 1);
    fprintf(stderr, "fft_x2k() %.3fs\n", (double)(clock()-t)/CLOCKS_PER_SEC);
}


void fft_k2x(fft_t *self, int offset_phase_off) {
    clock_t t = clock();
    fftw_execute(self->k2x);
    double L3inv = 1 / gsl_pow_3(self->L);  /* inverse box volume */
    for (long g = 0; g < self->Ng3_pad; ++g)
        self->f[g] *= L3inv;
    if (offset_phase_off)
        fft_offset_phase(self, -1);
    fprintf(stderr, "fft_k2x() %.3fs\n", (double)(clock()-t)/CLOCKS_PER_SEC);
}


void fft_offset_phase(fft_t *self, int onoff) {
    assert(onoff != 0);
    if (onoff > 0)
        onoff = 1;
    else
        onoff = -1;

    int Ng = self->Ng;
    double phase[3][Ng][2]; /* dim x mesh x complex */
    for (int i = 0; i < 3; ++i) {
        phase[i][0][0] = 1;
        phase[i][0][1] = 0;
    }
    for (int i = 0; i < 3; ++i)
    for (int j = 1; j <= Ng/2; ++j) {
        double arg = - onoff * 2*M_PI * j / Ng * self->offset[i];
        phase[i][j][0] = cos(arg);
        phase[i][j][1] = sin(arg);
        phase[i][Ng-j][0] = phase[i][j][0];
        phase[i][Ng-j][1] = - phase[i][j][1];
    }
    for (int i = 0; i < Ng; ++i)
    for (int j = 0; j < Ng; ++j)
    for (int k = 0; k <= Ng/2; ++k) {
        double phase_re = phase[0][i][0] * phase[1][j][0] * phase[2][k][0]
                        - phase[0][i][0] * phase[1][j][1] * phase[2][k][1]
                        - phase[0][i][1] * phase[1][j][0] * phase[2][k][1]
                        - phase[0][i][1] * phase[1][j][1] * phase[2][k][0];
        double phase_im = phase[0][i][1] * phase[1][j][0] * phase[2][k][0]
                        + phase[0][i][0] * phase[1][j][1] * phase[2][k][0]
                        + phase[0][i][0] * phase[1][j][0] * phase[2][k][1]
                        - phase[0][i][1] * phase[1][j][1] * phase[2][k][1];
        double field_re = F_Re(self,i,j,k);
        double field_im = F_Im(self,i,j,k);
        F_Re(self,i,j,k) = (field_re*phase_re - field_im*phase_im);
        F_Im(self,i,j,k) = (field_re*phase_im + field_im*phase_re);
    }
}


void fft_deconv(fft_t *self) {
    int Ng = self->Ng;
    double Winv[Ng];
    Winv[0] = 1;
    for (int i = 1; i <= Ng/2; ++i) {
        double arg = M_PI * i / Ng;
        Winv[i] = gsl_pow_4(arg / sin(arg));
        Winv[Ng-i] = Winv[i];
    }
    clock_t t = clock();
    for (int i = 0; i < Ng; ++i)
    for (int j = 0; j < Ng; ++j)
    for (int k = 0; k <= Ng/2; ++k) {
        double Winv3 = Winv[i] * Winv[j] * Winv[k];
        F_Re(self,i,j,k) *= Winv3;
        F_Im(self,i,j,k) *= Winv3;
    }
    fprintf(stderr, "fft_deconv() %.3fs\n", (double)(clock()-t)/CLOCKS_PER_SEC);
}


void fft_p2k(fft_t *self, gal_t *part) {
    fprintf(stderr, "fft_p2k() interlacing\n");

    double offset_dual[3] = {0.5, 0.5, 0.5};
    fft_p2g(self, part, offset_dual);
    fft_x2k(self, 1);

    double *f_dual = fft_exportf(self);

    double offset_none[3] = {0, 0, 0};
    fft_p2g(self, part, offset_none);
    fft_x2k(self, 0);

    for (long g = 0; g < self->Ng3_pad; ++g)
        self->f[g] = 0.5 * (f_dual[g] + self->f[g]);

    fft_deconv(self);

    free(f_dual);
}


double *fft_exportf(fft_t *self) {
    double *dest = (double *)malloc(self->Ng3_pad * sizeof(double)); assert(dest!=NULL);
    memcpy(dest, self->f, self->Ng3_pad * sizeof(double));
    return dest;
}


void fft_importf(fft_t *self, double *src) {
    assert(src != NULL);
    memcpy(self->f, src, self->Ng3_pad * sizeof(double));
}
