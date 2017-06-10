/* particles data type */

#ifndef GAL_H
#define GAL_H

#include <math.h>


/* squared distance with periodic boundary */
inline double pdist2(double x0, double y0, double z0,
                     double x1, double y1, double z1, double L){
    double dx = remainder(x1 - x0, L);
    double dy = remainder(y1 - y0, L);
    double dz = remainder(z1 - z0, L);
    return dx*dx + dy*dy + dz*dz;
}

/* normalize a vector */
inline void hat(double vec[3])
{
    double vecamp = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if(vecamp > 1e-7){
        vec[0] /= vecamp;
        vec[1] /= vecamp;
        vec[2] /= vecamp;
    }
}


typedef struct gal{
    long Np;
    double V;  /* volume */
    /* phase space */
    double *x, *y, *z;  /* [Mpc/h] */
    double *vx, *vy, *vz;  /* [km/s], optional */
    /* properties, optional */

    /* random catalog, optional */
    struct gal *rand;
} gal_t;

gal_t *gal_init(long Np, double V);

void gal_resize(gal_t *self, long Np);

void gal_free(gal_t *self);


/* move particles back into periodic boundary */
void gal_wrap(gal_t *self, double L);

/* move particles along LOS under global plane-parallel approximation
 * conformal Hubble a*H = a*H0*E(a), H0=100[h km/s/Mpc] */
gal_t *gal_rsd(gal_t *self, double los[3], double aH);


/* select Npsub particles inside a sub-box defined by
 * box = {x0, x1, y0, y1, z0, z1, L} */
gal_t *gal_subbox(gal_t *self, double box[7], double alpha);

/* select Npsub particles inside a sub-sphere defined by
 * sphere = {x0, y0, z0, R, L} */
gal_t *gal_subsphere(gal_t *self, double sphere[5], double alpha);


#endif
