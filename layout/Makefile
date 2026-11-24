CC ?= gcc
# -ffast-math lets float reductions vectorize; -march=native so both layouts
# get the same SIMD quality -- the comparison should be about memory, not codegen.
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -ffast-math -march=native

all: layout

layout: src/layout.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: layout
	./layout

clean:
	rm -f layout
