#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "../gal.h"
#include "io.h"


gal_t *gal_loadqpm_cubic(char *file, double L) {
    clock_t t = clock();
    FILE *fp = fopen(file, "r"); assert(fp!=NULL);
    int Np = 0, c;
    while ((c = fgetc(fp)) != EOF)
        if (c == '\n')
            ++ Np;
    rewind(fp);

    double V = L*L*L;
    gal_t *part = gal_init(Np, V);
    /* HACK */
    part->vx = (double *)malloc(Np * sizeof(double)); assert(part->vx!=NULL);
    part->vy = (double *)malloc(Np * sizeof(double)); assert(part->vy!=NULL);
    part->vz = (double *)malloc(Np * sizeof(double)); assert(part->vz!=NULL);
    for (int i = 0; i < Np; ++i) {
        int retval = fscanf(fp, "%lf %lf %lf %lf %lf %lf %*f %*d %*d\n",
                part->x+i, part->y+i, part->z+i,
                part->vx+i, part->vy+i, part->vz+i);
        assert(retval == 6);
    }
    fclose(fp);
    fprintf(stderr, "gal_loadqpm_cubic() %.2fs, loaded %d galaxies\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np);
    return part;
}
