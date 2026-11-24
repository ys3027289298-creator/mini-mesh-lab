CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: bake

bake: src/bake.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: bake
	./bake

clean:
	rm -f bake
