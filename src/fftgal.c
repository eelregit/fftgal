#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "fftgal.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


static double pow2(double x)
{
    return x*x;
}
static double pow3(double x)
{
    return x*x*x;
}


fftgal_t *fftgal_init(int Ng, double L, double V, int fold, char wisdom[])
{
    fftgal_t *self = (fftgal_t *)malloc(sizeof(fftgal_t)); assert(self!=NULL);
    self->Ng = Ng;
    self->L = L;
    if(V > 0)
        self->V = V;
    else
        self->V = pow3(L);
    self->fold = fold;
    self->Np3 = 0; /* until painted */
    for(int i=0; i<3; ++i)
        self->offset[i] = 0.; /* until painted */

    int Ng_pad = 2 * (Ng / 2 + 1);
    self->Ng_pad = Ng_pad;
    long int Ng3_pad = (long int)Ng * (long int)Ng * (long int)Ng_pad;
    self->Ng3_pad  = Ng3_pad;
    double *fx = (double *)fftw_malloc(Ng3_pad * sizeof(double)); assert(fx!=NULL);
    fftw_complex *fk = (fftw_complex *)fx;
    self->f = fx;

    if(!strcmp(wisdom, "FFTW_ESTIMATE")){
        fprintf(stderr, "fftgal_init() using FFTW_ESTIMATE\n");
        self->fx2fk = fftw_plan_dft_r2c_3d(Ng, Ng, Ng, fx, fk, FFTW_ESTIMATE);
        self->fk2fx = fftw_plan_dft_c2r_3d(Ng, Ng, Ng, fk, fx, FFTW_ESTIMATE);
    }
    else{
        int ret = fftw_import_wisdom_from_filename(wisdom);
        if(ret){
            fprintf(stderr, "fftgal_init() imported wisdom from %s\n", wisdom);
            self->fx2fk = fftw_plan_dft_r2c_3d(Ng, Ng, Ng, fx, fk, FFTW_ESTIMATE);
            self->fk2fx = fftw_plan_dft_c2r_3d(Ng, Ng, Ng, fk, fx, FFTW_ESTIMATE);
        }
        else{
            clock_t t = clock();
            self->fx2fk = fftw_plan_dft_r2c_3d(Ng, Ng, Ng, fx, fk, FFTW_MEASURE);
            self->fk2fx = fftw_plan_dft_c2r_3d(Ng, Ng, Ng, fk, fx, FFTW_MEASURE);
            ret = fftw_export_wisdom_to_filename(wisdom); assert(ret);
            fprintf(stderr, "fftgal_init() exported wisdom to %s\n", wisdom);
            fprintf(stderr, "fftgal_init() %.3fs to FFTW_MEASURE a %d^3 grid\n",
                    (double)(clock()-t)/CLOCKS_PER_SEC, Ng);
        }
    }

    return self;
}


void fftgal_x2fx(fftgal_t *self, double *x, double *y, double *z,
        long long int Np3, double offset[3])
{
    self->Np3 = Np3;
    int Ng = self->Ng;
    for(int i=0; i<3; ++i){
        offset[i] = remainder(offset[i], Ng);
        self->offset[i] = offset[i];
    }

    clock_t t = clock();
    for(long int g=0; g<self->Ng3_pad; ++g) /* must zero after fftw_plan */
        self->f[g] = 0.;
    double Hinv = Ng * self->fold / self->L;
    for(long long int p=0; p<Np3; ++p){
        double xp = x[p] * Hinv - offset[0];
        double yp = y[p] * Hinv - offset[1];
        double zp = z[p] * Hinv - offset[2];
        int xint = (int)floor(xp);
        int yint = (int)floor(yp);
        int zint = (int)floor(zp);
        double dx = xp - xint;
        double dy = yp - yint;
        double dz = zp - zint;
        xint = (xint % Ng + Ng) % Ng;
        yint = (yint % Ng + Ng) % Ng;
        zint = (zint % Ng + Ng) % Ng;
        int i[4] = {(xint - 1 + Ng) % Ng, xint, (xint + 1) % Ng, (xint + 2) % Ng};
        int j[4] = {(yint - 1 + Ng) % Ng, yint, (yint + 1) % Ng, (yint + 2) % Ng};
        int k[4] = {(zint - 1 + Ng) % Ng, zint, (zint + 1) % Ng, (zint + 2) % Ng};
        double wx[4] = {pow3(1 - dx) / 6, ((3*dx - 6)*pow2(dx) + 4) / 6,
                    (((-3*dx + 3)*dx + 3)*dx + 1) / 6, pow3(dx) / 6};
        double wy[4] = {pow3(1 - dy) / 6, ((3*dy - 6)*pow2(dy) + 4) / 6,
                    (((-3*dy + 3)*dy + 3)*dy + 1) / 6, pow3(dy) / 6};
        double wz[4] = {pow3(1 - dz) / 6, ((3*dz - 6)*pow2(dz) + 4) / 6,
                    (((-3*dz + 3)*dz + 3)*dz + 1) / 6, pow3(dz) / 6};
        for(int ii=0; ii<4; ++ii)
        for(int jj=0; jj<4; ++jj)
        for(int kk=0; kk<4; ++kk)
            F(self, i[ii], j[jj], k[kk]) += wx[ii] * wy[jj] * wz[kk];
    }
    double gpp = pow3(Ng * self->fold) / Np3; /* grid per particle, inverse density */
    gpp *= self->V / pow3(self->L);
    for(long int g=0; g<self->Ng3_pad; ++g)
        self->f[g] *= gpp;
    fprintf(stderr, "fftgal_x2fx() %.3fs to paint %lld particles to %d^3 grid\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3, Ng);
}


void fftgal_fx2fk(fftgal_t *self)
{
    clock_t t = clock();
    fftw_execute(self->fx2fk);
    double H3 = pow3(self->L / self->fold / self->Ng);
    for(long int g=0; g<self->Ng3_pad; ++g)
        self->f[g] *= H3;

    int is_offset = 0;
    for(int i=0; i<3; ++i)
        is_offset += fabs(self->offset[i]) > 1e-15;
    if(is_offset){
        int Ng = self->Ng;
        double phase[3][Ng][2]; /* dim x mesh x complex */
        for(int i=0; i<3; ++i){
            phase[i][0][0] = 1.;
            phase[i][0][1] = 0.;
        }
        for(int i=0; i<3; ++i)
        for(int j=1; j<=Ng/2; ++j){
            double arg = - 2 * M_PI * j / Ng * self->offset[i];
            phase[i][j][0] = cos(arg);
            phase[i][j][1] = sin(arg);
            phase[i][Ng-j][0] = phase[i][j][0];
            phase[i][Ng-j][1] = - phase[i][j][1];
        }
        for(int i=0; i<Ng; ++i)
        for(int j=0; j<Ng; ++j)
        for(int k=0; k<=Ng/2; ++k){
            double Rep = phase[0][i][0] * phase[1][j][0] * phase[2][k][0]
                       - phase[0][i][0] * phase[1][j][1] * phase[2][k][1]
                       - phase[0][i][1] * phase[1][j][0] * phase[2][k][1]
                       - phase[0][i][1] * phase[1][j][1] * phase[2][k][0];
            double Imp = phase[0][i][1] * phase[1][j][0] * phase[2][k][0]
                       + phase[0][i][0] * phase[1][j][1] * phase[2][k][0]
                       + phase[0][i][0] * phase[1][j][0] * phase[2][k][1]
                       - phase[0][i][1] * phase[1][j][1] * phase[2][k][1];
            double Red = F_Re(self,i,j,k);
            double Imd = F_Im(self,i,j,k);
            F_Re(self,i,j,k) = (Red*Rep - Imd*Imp);
            F_Im(self,i,j,k) = (Red*Imp + Imd*Rep);
        }
    }

    fprintf(stderr, "fftgal_fx2fk() %.3fs to FFT f(x) to f(k) on a %d^3 grid\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, self->Ng);
}


void fftgal_deconv(fftgal_t *self)
{
    int Ng = self->Ng;
    double *winv = (double *)malloc(sizeof(double) * Ng); assert(winv!=NULL);
    winv[0] = 1.;
    for(int i=1; i<=Ng/2; ++i){
        double arg = M_PI * i / Ng;
        winv[i] = pow2(pow2(arg / sin(arg)));
        winv[Ng-i] = winv[i];
    }
    clock_t t = clock();
    for(int i=0; i<Ng; ++i)
    for(int j=0; j<Ng; ++j)
    for(int k=0; k<=Ng/2; ++k){
        double winv3 = winv[i] * winv[j] * winv[k];
        F_Re(self,i,j,k) *= winv3;
        F_Im(self,i,j,k) *= winv3;
    }
    fprintf(stderr, "fftgal_deconv() %.3fs to deconvolve paintbrush on a %d^3 grid\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Ng);
    free(winv);
}


void fftgal_x2fk(fftgal_t *self, double *x, double *y, double *z, long long int Np3)
{
    fprintf(stderr, "fftgal_x2fk() interlacing with half-grid offset\n");
    double offset_dual[3] = {0.5, 0.5, 0.5};
    fftgal_x2fx(self, x, y, z, Np3, offset_dual);
    fftgal_fx2fk(self);

    double *fdual = fftgal_exportf(self);

    fprintf(stderr, "fftgal_x2fk() interlacing with zero offset\n");
    double offset_none[3] = {0., 0., 0.};
    fftgal_x2fx(self, x, y, z, Np3, offset_none);
    fftgal_fx2fk(self);

    for(long int g=0; g<self->Ng3_pad; ++g)
        self->f[g] = 0.5 * (fdual[g] + self->f[g]);
    fftgal_deconv(self);

    free(fdual);
}


void fftgal_fk2fx(fftgal_t *self)
{
    clock_t t = clock();
    fftw_execute(self->fk2fx);
    double Vinv = pow3(self->fold / self->L);
    for(long int g=0; g<self->Ng3_pad; ++g)
        self->f[g] *= Vinv;
    fprintf(stderr, "fftgal_fk2fx() %.3fs to FFT f(k) to f(x) on a %d^3 grid\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, self->Ng);
}


double *fftgal_exportf(fftgal_t *self)
{
    double *f_copy = (double *)malloc(self->Ng3_pad * sizeof(double));
    assert(f_copy!=NULL);
    memcpy(f_copy, self->f, self->Ng3_pad * sizeof(double));
    return f_copy;
}


void fftgal_importf(fftgal_t *self, double *f_copy)
{
    assert(f_copy!=NULL);
    memcpy(self->f, f_copy, self->Ng3_pad * sizeof(double));
}


void fftgal_savef(fftgal_t *self, char *filename)
{
    FILE *fp = fopen(filename, "wb"); assert(fp!=NULL);
    size_t ret = fwrite(self->f, sizeof(double), self->Ng3_pad, fp);
    assert(ret==self->Ng3_pad);
    fclose(fp);
}


void fftgal_readf(fftgal_t *self, char *filename)
{
    FILE *fp = fopen(filename, "rb"); assert(fp!=NULL);
    size_t ret = fread(self->f, sizeof(double), self->Ng3_pad, fp);
    assert(ret==self->Ng3_pad);
    fclose(fp);
}


void fftgal_free(fftgal_t *self)
{
    free(self->f);
    free(self);
}
