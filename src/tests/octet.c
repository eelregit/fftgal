#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "../fftgal.h"


int main()
{
    double eps = 1e-15;

    int Ng = 2;
    double L = 2.;
    fftgal_t *octet = fftgal_init(Ng, L, -1, 1, "FFTW_ESTIMATE");

    int Np = 1;
    double x[1] = {0.}, y[1] = {0.}, z[1] = {0.};
    double offset[3] = {7., 7., 7.};
    fftgal_x2fx(octet, x, y, z, Np, offset);
    fprintf(stdout, "* expected result: (2/3, 4/3)⨂ (2/3, 4/3)⨂ (2/3, 4/3)\n");
    for(int i=0; i<octet->Ng; ++i)
    for(int j=0; j<octet->Ng; ++j)
    for(int k=0; k<octet->Ng; ++k){
        double f = F(octet,i,j,k);
        double fexp = 8./27. * (i+1) * (j+1) * (k+1);
        double err = f - fexp;
        fprintf(stdout, "octet->fx[%d,%d,%d] = %.4f;  err = % .0e\n", i, j, k, f, err);
        assert(fabs(err) < eps);
    }

    fftgal_fx2fk(octet);
    fprintf(stdout, "* expected result: (2, 2/3)⨂ (2, 2/3)⨂ (2, 2/3)\n");
    for(int i=0; i<octet->Ng; ++i)
    for(int j=0; j<octet->Ng; ++j)
    for(int k=0; k<=octet->Ng/2; ++k){
        double f_re = F_Re(octet,i,j,k), f_im = F_Im(octet,i,j,k);
        double fexp_re = 8. / (1+2*i) / (1+2*j) / (1+2*k), fexp_im = 0.;
        double err_re = f_re - fexp_re, err_im = f_im - fexp_im;
        fprintf(stdout, "octet->fk[%d,%d,%d] = %.4f + %.4f i;", i, j, k, f_re, f_im);
        fprintf(stdout, "  err = % .0e + % .0e i\n", err_re, err_im);
        assert(fabs(err_re) < eps); assert(fabs(err_im) < eps);
    }

    fftgal_fk2fx(octet);
    fprintf(stdout, "* expected result: (4/3, 2/3)⨂ (4/3, 2/3)⨂ (4/3, 2/3)\n");
    for(int i=0; i<octet->Ng; ++i)
    for(int j=0; j<octet->Ng; ++j)
    for(int k=0; k<octet->Ng; ++k){
        double f = F(octet,i,j,k);
        double fexp = 64./27. / (i+1) / (j+1) / (k+1);
        double err = f - fexp;
        fprintf(stdout, "octet->fx[%d,%d,%d] = %.4f;  err = % .0e\n", i, j, k, f, err);
        assert(fabs(err) < eps);
    }

    fprintf(stdout, "* test PASSED on eps = %e\n", eps);
    fftgal_free(octet);
    return 0;
}
