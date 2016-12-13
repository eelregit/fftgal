#ifndef FFTGAL_H
#define FFTGAL_H

/* 3D cubic in-place two-way FFT from galaxies/particles painted with PCS,
 * providing interlacing interface, see Sefusatti et al. 2016
 */
#include <fftw3.h>


#define F(self,i,j,k) self->f[((long int)(k) + self->Ng_pad * ((long int)(j) + self->Ng * (long int)(i)))]


typedef struct {
    fftw_plan fx2fk;
    fftw_plan fk2fx;
    double *f; /* in-place: fk = (fftw_complex *)fx */
    int Ng;
    int Ng_pad; /* 2 * (Ng/2 + 1) */
    long int Ng3_pad; /* sizeof(f) / sizeof(double) */
    long long int Np3;
    double L;
} fftgal_t;


/* allocate f[] and plan fftw */
fftgal_t *fftgal_init(int Ng, double L, char wisdom[]);


/* paint , forward FFT, and deconvolve paintbrush */
void fftgal_x2fx(fftgal_t *self, double *x, double *y, double *z,
        long long int Np3, double offset);
void fftgal_fx2fk(fftgal_t *self);
void fftgal_deconv(fftgal_t *self);


/* simple interface to do interlaced painting+FFT, and average before deconvolve */
void fftgal_x2fk(fftgal_t *self, double *x, double *y, double *z, long long int Np3);


/* backward FFT */
void fftgal_fk2fx(fftgal_t *self);


/* make a copy of f[] */
double *fftgal_copyf(fftgal_t *self);


/* save and read f[] from files */
void fftgal_savef(fftgal_t *self, char *filename);
void fftgal_readf(fftgal_t *self, char *filename);


/* free fftgal */
void fftgal_kill(fftgal_t *self);

#endif
