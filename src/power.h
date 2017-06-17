#include "fft.h"
#include "gal.h"


/* P(k) = |delta(k)|^2 / V - Pshot
 * call this first */
void P(fft_t *grid, gal_t *part);

/* power spectrum multipoles
 * Pl's all proportional to P0 when los[] = {0,0,0} */
void Pl(fft_t *grid, gal_t *part, double los[3], double dK, char *file);
/* power spectrum wedges */
void Pmu(fft_t *grid, gal_t *part, double los[3], double dK, char *file);
