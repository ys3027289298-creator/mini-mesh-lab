CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: lightmap

lightmap: src/lightmap.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lightmap
	./lightmap

clean:
	rm -f lightmap
