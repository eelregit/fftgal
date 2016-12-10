CC = gcc
CFLAGS = -Wall

test_octet: test_octet.c fftgal.c fftgal.h
	$(CC) $(CFLAGS) test_octet.c fftgal.c -lfftw3 -lm -o test_octet
