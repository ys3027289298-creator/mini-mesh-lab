CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -pedantic
LDLIBS = -lm

.PHONY: all test export clean

all: uvseam

uvseam: src/uvseam.c
	$(CC) $(CFLAGS) -o $@ src/uvseam.c $(LDLIBS)

test: uvseam
	./uvseam --selftest

export: uvseam
	./uvseam --export out

clean:
	rm -f uvseam
	rm -rf out
