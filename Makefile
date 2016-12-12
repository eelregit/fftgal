CC = gcc
CFLAGS = -Wall

test_octet: test_octet.c fftgal.c fftgal.h
	$(CC) $(CFLAGS) test_octet.c fftgal.c -lfftw3 -lm -o test_octet

Plsb: Plsb.c fftgal.c fftgal.h power.c power.h box.c box.h io/qpm.c io/qpm.h
	$(CC) $(CFLAGS) Plsb.c fftgal.c power.c box.c io/qpm.c -lfftw3 -lm -o Plsb
