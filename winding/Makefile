CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: winding

winding: src/winding.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: winding
	./winding --out ../artifacts/winding_obj

clean:
	rm -f winding
	rm -rf roundtrip_out

.PHONY: all test clean
