#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


void qpm_cubic_mocks_read(char *filename, int N, double *x, double *y, double *z,
        double *vx, double *vy, double *vz, double *M, int *issat)
{
    x = (double *)malloc(N * sizeof(double)); assert(x!=NULL);
    y = (double *)malloc(N * sizeof(double)); assert(y!=NULL);
    z = (double *)malloc(N * sizeof(double)); assert(z!=NULL);
    vx = (double *)malloc(N * sizeof(double)); assert(vx!=NULL);
    vy = (double *)malloc(N * sizeof(double)); assert(vy!=NULL);
    vz = (double *)malloc(N * sizeof(double)); assert(vz!=NULL);
    M = (double *)malloc(N * sizeof(double)); assert(M!=NULL);
    issat = (int *)malloc(N * sizeof(int)); assert(issat!=NULL);

    FILE *fp = fopen(filename, "r"); assert(fp!=NULL);
    for(int i=0; i<N; ++i){
        int ret = fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf %d %*d\n",
                x+i, y+i, z+i, vx+i, vy+i, vz+i, M+i, issat+i);
        assert(ret==8);
    }
}
