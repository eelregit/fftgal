#include "fftgal.h"

/* under global plane-parallel approximation, move particles along some LOS
 * x, y, z, xd, yd, zd are in [Mpc/h]
 * vx, vy, vz are in [km/s]
 * conformal Hubble a*H = a*H0*E(a), H0=100[h km/s/Mpc]
 * if x, y, z are in [0, L]^3, xd, yd, zd may not
 * los[] doesn't have to be normalized
 * xd, yd, zd are allocated, remember to free them
 */
void rsd(double *x, double *y, double *z, double *vx, double *vy, double *vz,
        long long Np3, double **xd, double **yd, double **zd,
        double los[3], double aH);

/* multipoles, return number of bins
 * Pl all proportional to P0 when los[]={0,0,0} */
int Pl(fftgal_t *fg, double dK, double los[3], char *output);
