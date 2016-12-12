#include "fftgal.h"

/* under global plane-parallel approximation, move particles along some LOS
 * x, y, z, xd, yd, zd in [Mpc/h]
 * vx, vy, vz in [km/s]
 * the conformal Hubble a*H = a*H0*E(a), H0=100[h km/s/Mpc]
 * x, y, z in [0, L]^3 but xd, yd, zd may not
 * los[] doesn't have to be normalized
 */
void rsd(double *x, double *y, double *z, double *vx, double *vy, double *vz,
        long long int Np3, double *xd, double *yd, double *zd,
        double los[3], double aH);

/* multipoles, return number of bins
 * Pl all proportional to P0 when los[]={0,0,0} */
int Pl(fftgal_t *fg, double dK, double los[3], char *output);
