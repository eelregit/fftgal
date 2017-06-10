#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "../fft.h"
#include "../gal.h"


int main()
{
    double eps = 1e-15;

    int Ng = 2;
    double L = 2;
    fft_t *grid = fft_init(Ng, L, "FFTW_ESTIMATE");

    int Np = 1;
    double V = L*L*L;
    gal_t *part = gal_init(Np, V);
    part->x[0] = part->y[0] = part->z[0] = -8;

    double offset[3] = {7, 7, 7};
    fft_p2g(grid, part, offset);
    fprintf(stdout, "* expected result: (2/3, 4/3)⨂ (2/3, 4/3)⨂ (2/3, 4/3)\n");
    for(int i=0; i<Ng; ++i)
    for(int j=0; j<Ng; ++j)
    for(int k=0; k<Ng; ++k){
        double f = F(grid,i,j,k);
        double fexp = 8.0/27 * (i+1) * (j+1) * (k+1);
        double err = f - fexp;
        fprintf(stdout, "fx[%d,%d,%d] = %.4f;  err = % .0e\n", i, j, k, f, err);
        assert(fabs(err) < eps);
    }

    fft_x2k(grid, 1);
    fprintf(stdout, "* expected result: (2, 2/3)⨂ (2, 2/3)⨂ (2, 2/3)\n");
    for(int i=0; i<Ng; ++i)
    for(int j=0; j<Ng; ++j)
    for(int k=0; k<Ng/2; ++k){
        double f_re = F_Re(grid,i,j,k), f_im = F_Im(grid,i,j,k);
        double fexp_re = 8.0 / (1+2*i) / (1+2*j) / (1+2*k), fexp_im = 0;
        double err_re = f_re - fexp_re, err_im = f_im - fexp_im;
        fprintf(stdout, "fk[%d,%d,%d] = %.4f + %.4f i;", i, j, k, f_re, f_im);
        fprintf(stdout, "  err = % .0e + % .0e i\n", err_re, err_im);
        assert(fabs(err_re) < eps); assert(fabs(err_im) < eps);
    }

    fft_k2x(grid, 0);
    fprintf(stdout, "* expected result: (4/3, 2/3)⨂ (4/3, 2/3)⨂ (4/3, 2/3)\n");
    for(int i=0; i<Ng; ++i)
    for(int j=0; j<Ng; ++j)
    for(int k=0; k<Ng; ++k){
        double f = F(grid,i,j,k);
        double fexp = 64.0/27 / (i+1) / (j+1) / (k+1);
        double err = f - fexp;
        fprintf(stdout, "fx[%d,%d,%d] = %.4f;  err = % .0e\n", i, j, k, f, err);
        assert(fabs(err) < eps);
    }

    fprintf(stdout, "* test PASSED on eps = %e\n", eps);
    fft_free(grid);
    gal_free(part);
    return 0;
}
