/* move x, y, z back to [0, L) */
void pbc(double *x, double *y, double *z, long long int Np3, double L);

/* squared distance with periodic boundary */
double pdist2(double x1, double y1, double z1, double x2, double y2, double z2, double L);

/* select particles inside the subbox defined by
 * xyzlim={xmin, xmax, ymin, ymax, zmin, zmax}
 * Np3sb_max is a safe guess
 * return true Np3sb
 * xsb, ysb, zsb are allocated, remember to free them */
long long int subbox(double *x, double *y, double *z, long long int Np3, double xyzlim[6],
        double **xsb, double **ysb, double **zsb, long long int Np3sb_max);


/* select particles inside the sphere defined by
 * XYZRL={X, Y, Z, R, L}
 * fill the rest part of the bounding box with Poisson points
 * Np3ss_max is a safe guess
 * return Np3sb, rescaled from Np3ss so that the shot noise are kept the same
 * xsb, ysb, zsb are allocated, remember to free them */
long long int subsphere(double *x, double *y, double *z, long long int Np3, double XYZRL[5],
        double **xsb, double **ysb, double **zsb, long long int Np3ss_max);
