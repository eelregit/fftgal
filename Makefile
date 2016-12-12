CC = gcc
CFLAGS = -Wall -I$(FFTW_INC) -L$(FFTW_DIR)

.PHONY: clean

Plsb: Plsb.c fftgal.c fftgal.h power.c power.h box.c box.h io/qpm.c io/qpm.h
	$(CC) $(CFLAGS) Plsb.c fftgal.c power.c box.c io/qpm.c -lfftw3 -lm -o Plsb

test_octet: test_octet.c fftgal.c fftgal.h
	$(CC) $(CFLAGS) test_octet.c fftgal.c -lfftw3 -lm -o test_octet

clean:
	rm -f Plsb test_octet
