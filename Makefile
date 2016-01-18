CC = gcc
CFLAGS = -Wall -O3 -march=core2
LIBS = -lm
PREFIX = /usr/local/bin
SRC = src

ALL: swalign

clean:
	rm -f swalign 2>/dev/null

swalign: $(SRC)/swalign.c
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

install: swalign
	rm -f $(PREFIX)/swalign 2>/dev/null
	mv $^ $(PREFIX)/$^

remove:
	rm -f $(PREFIX)/swalign 2>/dev/null
