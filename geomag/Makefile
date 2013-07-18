CC = gcc
RM = rm -f
BIN = /usr/local/bin

all: magdec

magdec: geomag_calc.o magdec.o Makefile igrf11.h
	$(CC) -o $@ geomag_calc.o magdec.o -lm

install: magdec
	install magdec $(BIN)

clean:
	$(RM) *.o  magdec
