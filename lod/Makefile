CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: lod

lod: src/lod.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lod
	./lod

csv: lod
	./lod --csv

clean:
	rm -f lod
