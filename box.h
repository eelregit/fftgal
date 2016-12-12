#ifndef BOX_H
#define BOX_H

/* move x, y, z back to [0, L) */
void pbc(double *x, double *y, double *z, long long int Np3, double L);


/* save positions of particles inside the subbox defined by
 * xyzlim={xmin, xmax, ymin, ymax, zmin, zmax}
 * return Np3sb that must <= Np3sb_max */
long long int subbox(double *x, double *y, double *z, long long int Np3, double xyzlim[6],
        double *xsb, double *ysb, double *zsb, long long int Np3sb_max);

#endif
