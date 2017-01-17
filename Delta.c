/* measure subbox overdensity, it squared, and tide squared */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "fftgal.h"
#include "io/qpm.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


/* use bias≠1 to rescale to matter overdensity if it's known,
 * otherwise output is in tracer overdensity */
const double bias = 1.;


static double pow2(double x)
{
    return x*x;
}


int main(int argc, char *argv[])
{
    if(argc!=9){
        fprintf(stderr, "Usage: %s Ng L wisdom Nsb catdir a catid outdir\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    int Ng = atoi(argv[1]); assert(Ng>1 && Ng<=1024);
    double L = atof(argv[2]); assert(L>0. && L<1e4);
    char *wisdom = argv[3];
    int Nsb = atoi(argv[4]); assert(Nsb>0 && Nsb<=8);
    int Nsb3 = Nsb*Nsb*Nsb;
    char *catdir = argv[5];
    double a = atof(argv[6]); assert(a>0. && a<1.1);
    int catid = atoi(argv[7]); assert(catid>=1 && catid<=1000);
    char *outdir = argv[8];
    int Ngsb = Ng / Nsb; assert(Ng%Nsb==0);

    const int maxlen = 1024;
    char catalog[maxlen];
    int ret = snprintf(catalog, maxlen, "%s/a%.4f_%04d.mock", catdir, a, catid);
    assert(ret>=0 && ret<maxlen);
    double *x, *y, *z, *vx, *vy, *vz, *M;
    int *issat;
    int Np3 = qpm_cubic_mocks_read(catalog, &x, &y, &z, &vx, &vy, &vz, &M, &issat);
    assert(Np3>0);

    fftgal_t *fg = fftgal_init(Ng, L, wisdom);

    double *fk_copy = NULL;
    double Wamp[Ng], Wpha[3*Ng][2]; /* subbox smoothing */
    Wamp[0] = 1.;
    for(int i=1; i<Ng; ++i)
        Wamp[i] = sin(M_PI * i / Nsb) / sin(M_PI * i / Ng) / Ngsb;
    for(int i=0; i<3*Ng; ++i){
        Wpha[i][0] = cos(M_PI * i * (Ngsb-1) / Ng);
        Wpha[i][1] = sin(M_PI * i * (Ngsb-1) / Ng);
    }

    double Deltaij[9][Nsb3];
    for(int idim=0; idim<3; ++idim)
    for(int jdim=idim; jdim<3; ++jdim){
        if(fk_copy==NULL){
            fftgal_x2fx(fg, x, y, z, Np3, 0.);
            fftgal_fx2fk(fg);
            fk_copy = fftgal_exportf(fg);
        }
        else
            fftgal_importf(fg, fk_copy);

        fg->f[0] = 0.;
        for(int i=0; i<Ng; ++i){
            double Kvec[3];
            Kvec[0] = remainder(i, Ng);
            for(int j=0; j<Ng; ++j){
                Kvec[1] = remainder(j, Ng);
                for(int k=(i==0 && j==0); k<=Ng/2; ++k){ /* skip Kvec[]={0,0,0} */
                    Kvec[2] = k;
                    double Wamp3 = Wamp[i] * Wamp[j] * Wamp[k];
                    double ReW = Wpha[i+j+k][0] * Wamp3;
                    double ImW = Wpha[i+j+k][1] * Wamp3;
                    double Red = F(fg,i,j,2*k);
                    double Imd = F(fg,i,j,2*k+1);
                    double op = Kvec[idim] * Kvec[jdim]
                        / (Kvec[0]*Kvec[0] + Kvec[1]*Kvec[1] + Kvec[2]*Kvec[2]);
                    F(fg,i,j,2*k) = op * (Red*ReW - Imd*ImW) / bias;
                    F(fg,i,j,2*k+1) = op * (Red*ImW + Imd*ReW) / bias;
                }
            }
        }

        fftgal_fk2fx(fg);
        for(int isb=0; isb<Nsb; ++isb)
        for(int jsb=0; jsb<Nsb; ++jsb)
        for(int ksb=0; ksb<Nsb; ++ksb){
            Deltaij[3*idim+jdim][(isb*Nsb + jsb)*Nsb + ksb] = F(fg, isb*Ngsb, jsb*Ngsb, ksb*Ngsb);
            Deltaij[3*jdim+idim][(isb*Nsb + jsb)*Nsb + ksb] = F(fg, isb*Ngsb, jsb*Ngsb, ksb*Ngsb);
        }
    }
    double Delta[Nsb3], Delta2[Nsb3], S2[Nsb3];
    for(int isb=0; isb<Nsb; ++isb)
    for(int jsb=0; jsb<Nsb; ++jsb)
    for(int ksb=0; ksb<Nsb; ++ksb){
        int ijk_sb = (isb*Nsb + jsb)*Nsb + ksb;
        Delta[ijk_sb] = Deltaij[0][ijk_sb]
                                         + Deltaij[4][ijk_sb]
                                         + Deltaij[8][ijk_sb];
        Delta2[ijk_sb] = pow2(Delta2[ijk_sb]);
        S2[ijk_sb] = pow2(Deltaij[0][ijk_sb] - Delta[ijk_sb]/3.)
                   + pow2(Deltaij[1][ijk_sb]) + pow2(Deltaij[2][ijk_sb])
                   + pow2(Deltaij[3][ijk_sb]) + pow2(Deltaij[4][ijk_sb] - Delta[ijk_sb]/3.)
                   + pow2(Deltaij[5][ijk_sb]) + pow2(Deltaij[6][ijk_sb])
                   + pow2(Deltaij[7][ijk_sb]) + pow2(Deltaij[8][ijk_sb] - Delta[ijk_sb]/3.);
    }

    char outfile[maxlen];
    ret = snprintf(outfile, maxlen, "%s/a%.4f_%04d/Delta.txt", outdir, a, catid);
    assert(ret>=0 && ret<maxlen);
    FILE *fp = fopen(outfile, "w"); assert(fp!=NULL);
    fprintf(fp, "# ijk_sb Delta Delta^2 S_ij*S_ji\n");
    for(int isb=0; isb<Nsb; ++isb)
    for(int jsb=0; jsb<Nsb; ++jsb)
    for(int ksb=0; ksb<Nsb; ++ksb)
        fprintf(fp, "%d%d%d % e % e % e\n", isb, jsb, ksb,
                Delta[(isb*Nsb + jsb)*Nsb + ksb],
                Delta2[(isb*Nsb + jsb)*Nsb + ksb],
                S2[(isb*Nsb + jsb)*Nsb + ksb]);
    fclose(fp);

    free(x); free(y); free(z); free(vx); free(vy); free(vz); free(M); free(issat);
    free(fk_copy);
    fftgal_free(fg);
    return 0;
}
