# FFT-based analysis of large-scale structure of the Universe

* 3D, cubic, in-place FFT using FFTW
* paint with PCS scheme
* optional interlacing (Sefusatti et al. 2016)

* HACK: non-cubic geometry handled with self->V, and Poisson points
  outside V. Should be replaced with random catalog inside V
