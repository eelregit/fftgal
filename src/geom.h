/* move x, y, z back to [0, L) */
void pbc(double *x, double *y, double *z, long long Np3, double L);

/* squared distance with periodic boundary */
double pdist2(double x1, double y1, double z1, double x2, double y2, double z2, double L);

/* select Np3sub particles inside a sub-box defined by
 * xyzlim = {xmin, xmax, ymin, ymax, zmin, zmax}
 * fill the rest part of the bounding box with Np3rand=Np3bb-Np3sub Poisson points
 * bbox = {xbbmin, xbbmax, ybbmin, ybbmax, zbbmin, zbbmax} must contain the sub-region
 * pass bbox=xyzlim to turn off filling
 * Np3sub_max is a safe guess
 * return Np3bb, rescaled from Np3sub so that the shot noise are kept the same
 * xbb, ybb, zbb are allocated, remember to free them */
long long subbox(double *x, double *y, double *z, long long Np3, double xyzlim[6], double bbox[6],
        double **xbb, double **ybb, double **zbb, long long Np3sub_max);


/* select Np3sub particles inside a sub-sphere defined by
 * XYZRL = {X, Y, Z, R, L}
 * fill the rest part of the bounding box with Np3rand=Np3bb-Np3sub Poisson points
 * bbox = {xbbmin, xbbmax, ybbmin, ybbmax, zbbmin, zbbmax} must contain the sub-region
 * Np3sub_max is a safe guess
 * return Np3bb, rescaled from Np3sub so that the shot noise are kept the same
 * xbb, ybb, zbb are allocated, remember to free them */
long long subsphere(double *x, double *y, double *z, long long Np3, double XYZRL[5], double bbox[6],
        double **xbb, double **ybb, double **zbb, long long Np3sub_max);
