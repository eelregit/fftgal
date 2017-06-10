/* grid data type */

#ifndef FFT_H
#define FFT_H

#include <fftw3.h>
#include "gal.h"

/* field value in configuration space */
#define F(self,i,j,k) (self->f[((long)(k) + self->Ng_pad * ((long)(j) + self->Ng * (long)(i)))])
/* field value in Fourier space */
#define F_Re(self,i,j,k) F(self,i,j,2*k)
#define F_Im(self,i,j,k) F(self,i,j,2*k+1)


typedef struct fft{
    fftw_plan x2k, k2x;
    double *f;  /* in-place: fk = (fftw_complex *)fx */
    double offset[3];  /* grid (Dirac comb Ш) wrt boundary [0,L)^3 in grid unit */
    int Ng;
    int Ng_pad;  /* 2 * (Ng/2 + 1) */
    long Ng3_pad;  /* sizeof(f) / sizeof(double) */
    double L;  /* [Mpc/h] */
} fft_t;

fft_t *fft_init(int Ng, double L, char *wisdom);
void fft_free(fft_t *self);


void fft_p2g(fft_t *self, gal_t *part, double offset[3]);
double *fft_g2p(fft_t *self, gal_t *pos);


void fft_x2k(fft_t *self, int offset_phase_on);
void fft_k2x(fft_t *self, int offset_phase_off);
/* apply (onoff>0) or remove applied (onoff<0) offset phase */
void fft_offset_phase(fft_t *self, int onoff);


void fft_deconv(fft_t *self);

/* simple interlacing interface = (paint + FFT) x 2 + average + deconvolve */
void fft_p2k(fft_t *self, gal_t *part);


double *fft_exportf(fft_t *self);
void fft_importf(fft_t *self, double *src);


#endif
