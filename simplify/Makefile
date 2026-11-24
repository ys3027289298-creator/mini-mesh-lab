CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: simplify

simplify: src/simplify.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: simplify
	./simplify

clean:
	rm -f simplify
