CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: selfcheck

selfcheck: src/selfcheck.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: selfcheck
	./selfcheck

clean:
	rm -f selfcheck
