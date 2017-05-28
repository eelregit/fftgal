CC = gcc
CFLAGS = -Wall -I$(FFTW_INC) -L$(FFTW_DIR) -O2
LIBS = -lfftw3 -lm

.PHONY: clean

Plsb: src/app/Plsb.c src/fftgal.c src/fftgal.h src/power.c src/power.h \
    src/box.c src/box.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/Plsb.c src/fftgal.c src/power.c src/box.c \
		src/io/qpm.c $(LIBS) -o $@

Pl: src/app/Pl.c src/fftgal.c src/fftgal.h src/power.c src/power.h \
    src/box.c src/box.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/Pl.c src/fftgal.c src/power.c src/box.c \
		src/io/qpm.c $(LIBS) -o $@

ssm: src/app/ssm.c src/fftgal.c src/fftgal.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/ssm.c src/fftgal.c src/io/qpm.c $(LIBS) -o $@

Delta: src/app/Delta.c src/fftgal.c src/fftgal.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/Delta.c src/fftgal.c src/io/qpm.c $(LIBS) -o $@

test_octet: src/tests/octet.c src/fftgal.c src/fftgal.h
	$(CC) $(CFLAGS) src/tests/octet.c src/fftgal.c $(LIBS) -o $@

clean:
	rm -f Plsb Pl ssm Delta test_octet
