/* read Np3 galaxies, 8 columns skipping the halo id's
 * x, y, z, vx, vy, vz, M, issat are allocated, remember to free them */
void qpm_cubic_mocks_read(char *catalog, int Np3, double **x, double **y, double **z,
        double **vx, double **vy, double **vz, double **M, int **issat);
