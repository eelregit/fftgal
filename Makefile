CC = gcc
CFLAGS = -Wall -I$(FFTW_INC) -L$(FFTW_DIR) -O2

.PHONY: clean

Plsb: Plsb.c fftgal.c fftgal.h power.c power.h box.c box.h io/qpm.c io/qpm.h
	$(CC) $(CFLAGS) Plsb.c fftgal.c power.c box.c io/qpm.c -lfftw3 -lm -o Plsb

ssm: ssm.c fftgal.c fftgal.h io/qpm.c io/qpm.h
	$(CC) $(CFLAGS) ssm.c fftgal.c io/qpm.c -lfftw3 -lm -o ssm

Delta: Delta.c fftgal.c fftgal.h io/qpm.c io/qpm.h
	$(CC) $(CFLAGS) Delta.c fftgal.c io/qpm.c -lfftw3 -lm -o Delta

test_octet: test_octet.c fftgal.c fftgal.h
	$(CC) $(CFLAGS) test_octet.c fftgal.c -lfftw3 -lm -o test_octet

clean:
	rm -f Plsb ssm Delta test_octet
