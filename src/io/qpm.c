#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "../gal.h"
#include "io.h"


gal_t *gal_loadqpm_cubic(char *file, double L)
{
    FILE *fp = fopen(file, "r"); assert(fp!=NULL);
    int Np = 0, c;
    while((c = fgetc(fp)) != EOF)
        if(c=='\n')
            ++ Np;
    rewind(fp);

    gal_t *catalog = gal_init(Np, L*L*L);
    /* HACK */
    catalog->vx = (double *)malloc(Np * sizeof(double)); assert(catalog->vx!=NULL);
    catalog->vy = (double *)malloc(Np * sizeof(double)); assert(catalog->vy!=NULL);
    catalog->vz = (double *)malloc(Np * sizeof(double)); assert(catalog->vz!=NULL);

    clock_t t = clock();
    for(int i=0; i<Np; ++i){
        int retval = fscanf(fp, "%lf %lf %lf %lf %lf %lf %*f %*d %*d\n",
                catalog->x+i, catalog->y+i, catalog->z+i,
                catalog->vx+i, catalog->vy+i, catalog->vz+i);
        assert(retval==6);
    }
    fclose(fp);
    fprintf(stderr, "gal_loadqpm_cubic() %.3fs on loading %d galaxies\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np);

    return catalog;
}
