#include "fft.h"
#include "gal.h"


/* multipoles, return number of bins
 * Pl all proportional to P0 when los[] = {0,0,0} */
int Pl(fft_t *grid, gal_t *part, double dK, double los[3], char *file);
