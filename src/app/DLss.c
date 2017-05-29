/* measure super-sample modes from sub-spheres */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../fftgal.h"
#include "../io/qpm.h"


/* use bias!=1 to rescale to matter overdensity if it's known,
 * otherwise output is in tracer overdensity */
const double bias = 1.;


static double tophat(double KR) /* KR!=0 */
{
    return 3 * (sin(KR) - KR*cos(KR)) / gsl_pow_3(KR);
}


int main(int argc, char *argv[])
{
    if(argc!=9){
        fprintf(stderr, "Usage: %s Ng L wisdom Nss catdir a catid outdir\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    int Ng = atoi(argv[1]); assert(Ng>1 && Ng<=1024);
    double L = atof(argv[2]); assert(L>0. && L<1e4);
    char *wisdom = argv[3];
    int Nss = atoi(argv[4]); assert(Nss>0 && Nss<=8);
    char *catdir = argv[5];
    double a = atof(argv[6]); assert(a>0. && a<1.1);
    int catid = atoi(argv[7]); assert(catid>=1 && catid<=1000);
    char *outdir = argv[8];
    int Ngss = Ng / Nss; assert(Ng%Nss==0);

    const int maxlen = 1024;
    char catalog[maxlen];
    int ret = snprintf(catalog, maxlen, "%s/a%.4f_%04d.mock", catdir, a, catid);
    assert(ret>=0 && ret<maxlen);
    double *x, *y, *z, *vx, *vy, *vz, *M;
    int *issat;
    int Np3 = qpm_cubic_mocks_read(catalog, &x, &y, &z, &vx, &vy, &vz, &M, &issat);
    assert(Np3>0);

    fftgal_t *fg = fftgal_init(Ng, L, -1, 1, wisdom);

    double *fk_copy = NULL;
    for(int ilos=0; ilos<2; ++ilos)
    for(int jlos=0; jlos<2; ++jlos)
    for(int klos=0+(ilos==0 && jlos==0); klos<2-(ilos==1 && jlos==1); ++klos){ /* skip {0,0,0}, {1,1,1} */
        fprintf(stderr, "================== los div ==================\n");
        double losamp = sqrt(ilos*ilos + jlos*jlos + klos*klos);
        double loshat[3];
        loshat[0] = ilos / losamp;
        loshat[1] = jlos / losamp;
        loshat[2] = klos / losamp;

        double DeltaL[3*Nss*Nss*Nss];
        for(int L_half=0; L_half<=2; ++L_half){
            if(fk_copy==NULL){
                double offset[3] = {0., 0., 0.};
                fftgal_x2fx(fg, x, y, z, Np3, offset);
                fftgal_fx2fk(fg);
                fk_copy = fftgal_exportf(fg);
            }
            else
                fftgal_importf(fg, fk_copy);

            double Kval[Ng];
            Kval[0] = 0.;
            for(int i=1; i<=Ng/2; ++i){
                Kval[i] = i;
                Kval[Ng-i] = - i;
            }
            fg->f[0] = 0.;
            for(int i=0; i<Ng; ++i)
            for(int j=0; j<Ng; ++j)
            for(int k=(i==0 && j==0); k<=Ng/2; ++k){ /* skip {0,0,0} */
                double Kamp = sqrt(gsl_pow_2(Kval[i]) + gsl_pow_2(Kval[j]) + gsl_pow_2(Kval[k]));
                double mu2 = gsl_pow_2((Kval[i]*loshat[0] + Kval[j]*loshat[1] + Kval[k]*loshat[2]) / Kamp);
                double Legendre[3] = {1., 1.5*mu2 - 0.5, (4.375*mu2 - 3.75)*mu2 + 0.375};
                double W = tophat(Kamp * 2*M_PI / Nss);
                F_Re(fg,i,j,k) *= Legendre[L_half] * W;
                F_Im(fg,i,j,k) *= Legendre[L_half] * W;
            }

            fftgal_fk2fx(fg);
            for(int iss=0; iss<Nss; ++iss)
            for(int jss=0; jss<Nss; ++jss)
            for(int kss=0; kss<Nss; ++kss)
                DeltaL[((L_half*Nss + iss)*Nss + jss)*Nss + kss] = F(fg, iss*Ngss, jss*Ngss, kss*Ngss) / bias;
        }

        char outfile[maxlen];
        ret = snprintf(outfile, maxlen, "%s/a%.4f_%04d/DLss_los%d%d%d.txt",
                outdir, a, catid, ilos, jlos, klos);
        assert(ret>=0 && ret<maxlen);
        FILE *fp = fopen(outfile, "w"); assert(fp!=NULL);
        fprintf(fp, "# ijk_ss Delta_0 Delta_2 Delta_4\n");
        for(int iss=0; iss<Nss; ++iss)
        for(int jss=0; jss<Nss; ++jss)
        for(int kss=0; kss<Nss; ++kss)
            fprintf(fp, "%d%d%d % e % e % e\n", iss, jss, kss,
                    DeltaL[((0*Nss + iss)*Nss + jss)*Nss + kss],
                    DeltaL[((1*Nss + iss)*Nss + jss)*Nss + kss],
                    DeltaL[((2*Nss + iss)*Nss + jss)*Nss + kss]);
        fclose(fp);
    }

    free(x); free(y); free(z); free(vx); free(vy); free(vz); free(M); free(issat);
    free(fk_copy);
    fftgal_free(fg);
    return 0;
}
