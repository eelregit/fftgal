CC = gcc
CFLAGS = -Wall -I$(FFTW_INC) -L$(FFTW_DIR) -O2
LIBS = -lfftw3 -lm

.PHONY: clean

Pl: src/app/Pl.c src/fftgal.c src/fftgal.h src/power.c src/power.h \
    src/geom.c src/geom.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/Pl.c src/fftgal.c src/power.c src/geom.c \
		src/io/qpm.c $(LIBS) -o $@

Plsb: src/app/Plsb.c src/fftgal.c src/fftgal.h src/power.c src/power.h \
    src/geom.c src/geom.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/Plsb.c src/fftgal.c src/power.c src/geom.c \
		src/io/qpm.c $(LIBS) -o $@

Plss: src/app/Plss.c src/fftgal.c src/fftgal.h src/power.c src/power.h \
    src/geom.c src/geom.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/Plss.c src/fftgal.c src/power.c src/geom.c \
		src/io/qpm.c $(LIBS) -o $@

DLsb: src/app/DLsb.c src/fftgal.c src/fftgal.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/DLsb.c src/fftgal.c src/io/qpm.c $(LIBS) -o $@

DLss: src/app/DLss.c src/fftgal.c src/fftgal.h src/io/qpm.c src/io/qpm.h
	$(CC) $(CFLAGS) src/app/DLss.c src/fftgal.c src/io/qpm.c $(LIBS) -o $@

test_octet: src/tests/octet.c src/fftgal.c src/fftgal.h
	$(CC) $(CFLAGS) src/tests/octet.c src/fftgal.c $(LIBS) -o $@

clean:
	rm -f Pl Plsb Plss DLsb DLss Delta test_octet
