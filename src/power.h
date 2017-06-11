#include "fft.h"
#include "gal.h"


/* power spectrum multipoles
 * Pl's all proportional to P0 when los[] = {0,0,0} */
void Pl(fft_t *grid, gal_t *part, double los[3], double dK, char *file);
