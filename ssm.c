/* measure super-sample modes (mean and tide) from subboxes */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "fftgal.h"
#include "io/qpm.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


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
    for(int ilos=0; ilos<2; ++ilos)
    for(int jlos=0; jlos<2; ++jlos)
    for(int klos=0; klos<2-(ilos==1 && jlos==1); ++klos){ /* skip los[] = {1,1,1} */
        fprintf(stderr, "================== los div ==================\n");
        double losamp = sqrt(ilos*ilos + jlos*jlos + klos*klos);
        double loshat[3];
        if(losamp > 1e-7){
            loshat[0] = ilos / losamp;
            loshat[1] = jlos / losamp;
            loshat[2] = klos / losamp;
        }
        else
            loshat[0] = loshat[1] = loshat[2] = 0.;

        if(fk_copy==NULL){
            fftgal_x2fx(fg, x, y, z, Np3, 0.);
            fftgal_fx2fk(fg);
            fk_copy = fftgal_exportf(fg);
        }
        else
            fftgal_importf(fg, fk_copy);

        double Wamp[Ng], Wpha[3*Ng][2]; /* subbox smoothing */
        Wamp[0] = 1.;
        for(int i=1; i<Ng; ++i)
            Wamp[i] = sin(M_PI * i / Nsb) / sin(M_PI * i / Ng) / Ngsb;
        for(int i=0; i<3*Ng; ++i){
            Wpha[i][0] = cos(M_PI * i * (Ngsb-1) / Ng);
            Wpha[i][1] = sin(M_PI * i * (Ngsb-1) / Ng);
        }
        fg->f[0] = 0.;
        for(int i=0; i<Ng; ++i){
            double Kvec[3];
            Kvec[0] = remainder(i, Ng);
            for(int j=0; j<Ng; ++j){
                Kvec[1] = remainder(j, Ng);
                for(int k=(i==0 && j==0); k<=Ng/2; ++k){ /* skip Kvec[]={0,0,0} */
                    Kvec[2] = k;
                    double mu2 = pow2(Kvec[0]*loshat[0] + Kvec[1]*loshat[1] + Kvec[2]*loshat[2])
                            / (Kvec[0]*Kvec[0] + Kvec[1]*Kvec[1] + Kvec[2]*Kvec[2]);
                    double Legendre = ilos+jlos+klos==0 ? 1. : 1.5*mu2-0.5;
                    double Wamp3 = Wamp[i] * Wamp[j] * Wamp[k];
                    double ReW = Wpha[i+j+k][0] * Wamp3;
                    double ImW = Wpha[i+j+k][1] * Wamp3;
                    double Red = F(fg,i,j,2*k);
                    double Imd = F(fg,i,j,2*k+1);
                    F(fg,i,j,2*k) = Legendre * (Red*ReW - Imd*ImW);
                    F(fg,i,j,2*k+1) = Legendre * (Red*ImW + Imd*ReW);
                }
            }
        }

        fftgal_fk2fx(fg);

        char outfile[maxlen];
        ret = snprintf(outfile, maxlen, "%s/a%.4f_%04d/ssm_los%d%d%d.txt",
                outdir, a, catid, ilos, jlos, klos);
        assert(ret>=0 && ret<maxlen);
        FILE *fp = fopen(outfile, "w"); assert(fp!=NULL);
        if(ilos+jlos+klos == 0)
            fprintf(fp, "# isb jsb ksb Delta_0\n");
        else
            fprintf(fp, "# isb jsb ksb Delta_2\n");
        for(int isb=0; isb<Nsb; ++isb)
        for(int jsb=0; jsb<Nsb; ++jsb)
        for(int ksb=0; ksb<Nsb; ++ksb)
            fprintf(fp, "%d %d %d % e\n", isb, jsb, ksb, F(fg, isb*Ngsb, jsb*Ngsb, ksb*Ngsb));
        fclose(fp);
    }

    free(x); free(y); free(z); free(vx); free(vy); free(vz); free(M); free(issat);
    free(fk_copy);
    fftgal_free(fg);
    return 0;
}
