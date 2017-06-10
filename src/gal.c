#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "gal.h"


gal_t *gal_init(long Np, double V){
    gal_t *self = (gal_t *)malloc(sizeof(gal_t)); assert(self!=NULL);
    self->Np = Np;
    self->V = V;
    self->x = (double *)malloc(self->Np * sizeof(double)); assert(self->x!=NULL);
    self->y = (double *)malloc(self->Np * sizeof(double)); assert(self->y!=NULL);
    self->z = (double *)malloc(self->Np * sizeof(double)); assert(self->z!=NULL);
    self->vx = self->vy = self->vz = NULL;
    self->rand = NULL;
    return self;
}


void gal_resize(gal_t *self, long Np){
    self->Np = Np;
    self->x = (double *)realloc(self->x, Np * sizeof(double)); assert(self->x!=NULL);
    self->y = (double *)realloc(self->y, Np * sizeof(double)); assert(self->y!=NULL);
    self->z = (double *)realloc(self->z, Np * sizeof(double)); assert(self->z!=NULL);
    if(self->vx != NULL){
        self->vx = (double *)realloc(self->vx, Np * sizeof(double)); assert(self->vx!=NULL);
    }
    if(self->vy != NULL){
        self->vy = (double *)realloc(self->vy, Np * sizeof(double)); assert(self->vy!=NULL);
    }
    if(self->vz != NULL){
        self->vz = (double *)realloc(self->vz, Np * sizeof(double)); assert(self->vz!=NULL);
    }
}


void gal_free(gal_t *self){
    if(self != NULL){
        gal_free(self->rand);
        free(self->x); free(self->y); free(self->z);
        free(self->vx); free(self->vy); free(self->vz);
    }
    free(self);
}


void gal_wrap(gal_t *self, double L)
{
    double Linv = 1 / L;
    clock_t t = clock();
    for(long p=0; p<self->Np; ++p){
        self->x[p] -= L * floor(self->x[p] * Linv);
        self->y[p] -= L * floor(self->y[p] * Linv);
        self->z[p] -= L * floor(self->z[p] * Linv);
    }
    fprintf(stderr, "gal_wrap() %.3fs\n", (double)(clock()-t)/CLOCKS_PER_SEC);
}


gal_t *gal_rsd(gal_t *self, double los[3], double aH)
{
    gal_t *distorted = gal_init(self->Np, self->V);
    hat(los);
    double aHinv = 1 / aH;
    clock_t t = clock();
    for(long p=0; p<distorted->Np; ++p){
        double vdotlos = los[0]*self->vx[p] + los[1]*self->vy[p] + los[2]*self->vz[p];
        distorted->x[p] = self->x[p] + aHinv * vdotlos * los[0];
        distorted->y[p] = self->y[p] + aHinv * vdotlos * los[1];
        distorted->z[p] = self->z[p] + aHinv * vdotlos * los[2];
    }
    fprintf(stderr, "gal_rsd() %.3fs, los[]={%.3f,%.3f,%.3f}\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, los[0], los[1], los[2]);
    return distorted;
}


gal_t *gal_subbox(gal_t *self, double box[7], double alpha)
{
    double x0 = box[0], x1 = box[1];
    double y0 = box[2], y1 = box[3];
    double z0 = box[4], z1 = box[5];
    double L = box[6];
    x1 = x0 + remainder(x1-x0, L); assert(x1>x0);
    y1 = y0 + remainder(y1-y0, L); assert(y1>y0);
    z1 = z0 + remainder(z1-z0, L); assert(z1>z0);
    double Vsub = (x1-x0) * (y1-y0) * (z1-z0);

    gal_t *sub = gal_init(self->Np, Vsub);
    clock_t t = clock();
    long p, Npsub;
    for(p=0, Npsub=0; p<self->Np; ++p)
        if(remainder(self->x[p] - x0, L)>=0 && remainder(self->x[p] - x1, L)<0
        && remainder(self->y[p] - y0, L)>=0 && remainder(self->y[p] - y1, L)<0
        && remainder(self->z[p] - z0, L)>=0 && remainder(self->z[p] - z1, L)<0){
            sub->x[Npsub] = self->x[p];
            sub->y[Npsub] = self->y[p];
            sub->z[Npsub] = self->z[p];
            ++ Npsub;
        }
    gal_resize(sub, Npsub);

    long Nprand = lround(Npsub / alpha);
    gal_t *rand = gal_init(Nprand, Vsub);
    static gsl_rng *rng = NULL;
    if(rng == NULL){
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_default);
    }
    for(p=0; p<Nprand; ++p){
        rand->x[p] = x0 + (x1-x0) * gsl_rng_uniform(rng);
        rand->y[p] = y0 + (y1-y0) * gsl_rng_uniform(rng);
        rand->z[p] = z0 + (z1-z0) * gsl_rng_uniform(rng);
    }
    sub->rand = rand;

    fprintf(stderr, "gal_subbox() %.3fs on picking out %ld+%ld particles\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Npsub, Nprand);
    return sub;
}


gal_t *gal_subsphere(gal_t *self, double sphere[5], double alpha)
{
    double x0 = sphere[0], y0 = sphere[1], z0 = sphere[2];
    double R = sphere[3], L = sphere[4];
    double R2 = R*R, Vsub = 4*M_PI/3 * R2*R;
    assert(2*R <= L);

    gal_t *sub = gal_init(self->Np, Vsub);
    clock_t t = clock();
    long p, Npsub;
    for(p=0, Npsub=0; p<self->Np; ++p)
        if(pdist2(self->x[p], self->y[p], self->z[p], x0, y0, z0, L) < R2){
            sub->x[Npsub] = self->x[p];
            sub->y[Npsub] = self->y[p];
            sub->z[Npsub] = self->z[p];
            ++ Npsub;
        }
    gal_resize(sub, Npsub);

    long Nprand = lround(Npsub / alpha);
    gal_t *rand = gal_init(Nprand, Vsub);
    static gsl_rng *rng = NULL;
    if(rng == NULL){
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_default);
    }
    for(p=0; p<Nprand; ++p){
        do{
            rand->x[p] = x0 + R * (2*gsl_rng_uniform(rng) - 1);
            rand->y[p] = y0 + R * (2*gsl_rng_uniform(rng) - 1);
            rand->z[p] = z0 + R * (2*gsl_rng_uniform(rng) - 1);
        }while(pdist2(rand->x[p], rand->y[p], rand->z[p], x0, y0, z0, L) >= R2);
    }
    sub->rand = rand;

    fprintf(stderr, "gal_subsphere() %.3fs on picking out %ld+%ld particles\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Npsub, Nprand);
    return sub;
}
