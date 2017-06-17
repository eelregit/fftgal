CC = gcc
CFLAGS = -Wall -I$(FFTW_INC) -L$(FFTW_DIR) $(shell gsl-config --cflags) -O2
LIBS = -lfftw3 $(shell gsl-config --libs)

.PHONY: all clean

all: test_octet P Psb Pss DLsb DLss

clean:
	rm -f test_octet P Psb Pss DLsb DLss

test_octet: src/tests/octet.c src/fft.c src/fft.h src/gal.c src/gal.h
	$(CC) $(CFLAGS) src/tests/octet.c src/fft.c src/gal.c $(LIBS) -o $@

P: src/app/P.c src/fft.c src/fft.h src/gal.c src/gal.h \
	src/power.c src/power.h src/io/qpm.c src/io/io.h
	$(CC) $(CFLAGS) src/app/P.c src/fft.c src/gal.c src/power.c \
		src/io/qpm.c $(LIBS) -o $@

Psb: src/app/Psb.c src/fft.c src/fft.h src/gal.c src/gal.h \
	src/power.c src/power.h src/io/qpm.c src/io/io.h
	$(CC) $(CFLAGS) src/app/Psb.c src/fft.c src/gal.c src/power.c \
		src/io/qpm.c $(LIBS) -o $@

Pss: src/app/Pss.c src/fft.c src/fft.h src/gal.c src/gal.h \
	src/power.c src/power.h src/io/qpm.c src/io/io.h
	$(CC) $(CFLAGS) src/app/Pss.c src/fft.c src/gal.c src/power.c \
		src/io/qpm.c $(LIBS) $(GSL_FLAGS) -o $@

DLsb: src/app/DLsb.c src/fft.c src/fft.h src/gal.c src/gal.h \
	src/io/qpm.c src/io/io.h
	$(CC) $(CFLAGS) src/app/DLsb.c src/fft.c src/gal.c src/io/qpm.c $(LIBS) -o $@

DLss: src/app/DLss.c src/fft.c src/fft.h src/gal.c src/gal.h \
	src/io/qpm.c src/io/io.h
	$(CC) $(CFLAGS) src/app/DLss.c src/fft.c src/gal.c src/io/qpm.c $(LIBS) -o $@
