/* under global plane-parallel approximation, move particles along some LOS
 * x, y, z, xd, yd, zd in [Mpc/h]
 * vx, vy, vz in [km/s]
 * the conformal Hubble a*H = a*H0*E(a), H0=100[h km/s/Mpc]
 * x, y, z in [0, L]^3 but xd, yd, zd may not
 */
void rsd(double *x, double *y, double *z, double *vx, double *vy, double *vz,
        long long int Np3, double *xd, double *yd, double *zd,
        double los[3], double aH);
