#include <math.h>
#include <fenv.h>


void pbc(double *x, double *y, double *z, long long int Np3, double L)
{
    double Linv = 1 / L;
    int ret = fesetround(FE_DOWNWARD); assert(!ret);
    long long int p;
    for(p=0; p<Np3; ++p){
        x[p] -= L * nearbyint(x[p] * Linv);
        y[p] -= L * nearbyint(y[p] * Linv);
        z[p] -= L * nearbyint(z[p] * Linv);
    }
}
