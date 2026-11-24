CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: cull

cull: src/cull.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: cull
	./cull

clean:
	rm -f cull
