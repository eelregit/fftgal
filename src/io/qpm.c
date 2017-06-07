#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "qpm.h"


int qpm_cubic_mocks_load(char *catalog, double **x, double **y, double **z,
        double **vx, double **vy, double **vz, double **M, int **issat)
{
    FILE *fp = fopen(catalog, "r"); assert(fp!=NULL);
    int Np3 = 0, ch;
    while((ch=fgetc(fp)) != EOF)
        if(ch=='\n')
            ++ Np3;
    rewind(fp);
    *x = (double *)malloc(Np3 * sizeof(double)); assert(*x!=NULL);
    *y = (double *)malloc(Np3 * sizeof(double)); assert(*y!=NULL);
    *z = (double *)malloc(Np3 * sizeof(double)); assert(*z!=NULL);
    *vx = (double *)malloc(Np3 * sizeof(double)); assert(*vx!=NULL);
    *vy = (double *)malloc(Np3 * sizeof(double)); assert(*vy!=NULL);
    *vz = (double *)malloc(Np3 * sizeof(double)); assert(*vz!=NULL);
    *M = (double *)malloc(Np3 * sizeof(double)); assert(*M!=NULL);
    *issat = (int *)malloc(Np3 * sizeof(int)); assert(*issat!=NULL);
    clock_t t = clock();
    for(int i=0; i<Np3; ++i){
        int ret = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %d %*d\n",
                *x+i, *y+i, *z+i, *vx+i, *vy+i, *vz+i, *M+i, *issat+i);
        assert(ret==8);
    }
    fclose(fp);
    fprintf(stderr, "qpm_cubic_mocks_load() %.3fs to load %d galaxies\n",
            (double)(clock()-t)/CLOCKS_PER_SEC, Np3);
    return Np3;
}
