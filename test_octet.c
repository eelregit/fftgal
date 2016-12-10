#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "fftgal.h"

int main()
{
    fftgal_t octet = fftgal_init(2, 2., "FFTW_ESTIMATE");

    double x[1] = {0.}, y[1] = {0.}, z[1] = {0.};
    fftgal_x2fx(octet, x, y, z, 1, 0.5);

    int i, j, k;
    for(i=0; i<octet.Ng; ++i)
        for(j=0; j<octet.Ng; ++j)
            for(k=0; k<octet.Ng; ++k)
                fprintf(stdout, "octet.fx[%d,%d,%d] = %.16f\n", i, j, k, F(octet,i,j,k));
    double err = F(octet,0,0,0)-64./27;
    double eps = 1e-14;
    fprintf(stdout, "* error on octet.fx[0,0,0] = %e\n", err); assert(fabs(err) < eps);
    /* expected result: 8/27*{8, 4, 4, 4, 2, 2, 2, 1} */

    fftgal_fx2fk(octet);
    for(i=0; i<octet.Ng; ++i)
        for(j=0; j<octet.Ng; ++j)
            for(k=0; k<=octet.Ng/2; ++k)
                fprintf(stdout, "octet.fk[%d,%d,%d] = %.16f + %.16f i\n", i, j, k,
                        F(octet,i,j,2*k), F(octet,i,j,2*k+1));
    err = F(octet,1,1,2)-8./27;
    fprintf(stdout, "* error on Re(octet.fk[1,1,1]) = %e\n", err); assert(fabs(err) < eps);
    err = F(octet,1,1,3)-0;
    fprintf(stdout, "* error on Im(octet.fk[1,1,1]) = %e\n", err); assert(fabs(err) < eps);

    fftgal_fk2fx(octet);
    for(i=0; i<octet.Ng; ++i)
        for(j=0; j<octet.Ng; ++j)
            for(k=0; k<octet.Ng; ++k)
                fprintf(stdout, "octet.fx[%d,%d,%d] = %.16f\n", i, j, k, F(octet,i,j,k));
    err = F(octet,1,1,0)-16./27;
    fprintf(stdout, "* error on octet.fx[1,1,0] = %e\n", err); assert(fabs(err) < eps);

    fprintf(stdout, "* test_octet PASSED on eps = %e\n", eps);
    fftgal_freef(octet);
    return 0;
}
