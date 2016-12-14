/* read Np3 galaxies, 8 columns skipping the halo id's
 * x, y, z, vx, vy, vz, M, issat are allocated, remember to free them */
int qpm_cubic_mocks_read(char *catalog, double **x, double **y, double **z,
        double **vx, double **vy, double **vz, double **M, int **issat);
