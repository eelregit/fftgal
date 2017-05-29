#ifndef FFTGAL_H
#define FFTGAL_H

#include <fftw3.h>


#define F(self,i,j,k) (self->f[((long int)(k) + self->Ng_pad * ((long int)(j) + self->Ng * (long int)(i)))])
#define F_Re(self,i,j,k) F(self,i,j,2*k)
#define F_Im(self,i,j,k) F(self,i,j,2*k+1)


typedef struct {
    fftw_plan fx2fk;
    fftw_plan fk2fx;
    double *f; /* in-place: fk = (fftw_complex *)fx */
    int Ng;
    int Ng_pad; /* 2 * (Ng/2 + 1) */
    long int Ng3_pad; /* sizeof(f) / sizeof(double) */
    long long int Np3;
    double L;
    double V; /* L^3 or volume for non-cubic geometry */
    int fold; /* folded length is L/fold */
    double offset[3]; /* mesh (Dirac comb Ш) wrt boundary */
} fftgal_t;


/* allocate f[] and plan fftw
 * set V<0 for cubes
 * set 0<V<L^3 for other geometry bounded by the cube */
fftgal_t *fftgal_init(int Ng, double L, double V, int fold, char wisdom[]);


/* paint , forward FFT, and deconvolve paintbrush */
void fftgal_x2fx(fftgal_t *self, double *x, double *y, double *z,
        long long int Np3, double offset[3]);
void fftgal_fx2fk(fftgal_t *self);
void fftgal_deconv(fftgal_t *self);


/* simple interface to do interlaced painting+FFT, and average before deconvolve */
void fftgal_x2fk(fftgal_t *self, double *x, double *y, double *z, long long int Np3);


/* backward FFT */
void fftgal_fk2fx(fftgal_t *self);


/* export and import a copy of f[] */
double *fftgal_exportf(fftgal_t *self);
void fftgal_importf(fftgal_t *self, double *f_copy);


/* save and read f[] from files */
void fftgal_savef(fftgal_t *self, char *filename);
void fftgal_readf(fftgal_t *self, char *filename);


/* free fftgal */
void fftgal_free(fftgal_t *self);

#endif
