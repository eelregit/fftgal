/* move x, y, z back to [0, L) */
void pbc(double *x, double *y, double *z, long long int Np3, double L);


/* select particles inside the subbox defined by
 * xyzlim={xmin, xmax, ymin, ymax, zmin, zmax}
 * Np3sb_max is a safe guess
 * return true Np3sb
 * xsb, ysb, zsb are allocated, remember to free them */
long long int subbox(double *x, double *y, double *z, long long int Np3, double xyzlim[6],
        double **xsb, double **ysb, double **zsb, long long int Np3sb_max);
