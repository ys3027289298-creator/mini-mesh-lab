CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: triangulate

triangulate: src/triangulate.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: triangulate
	./triangulate

clean:
	rm -f triangulate
