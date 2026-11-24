CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: holefill

holefill: src/holefill.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: holefill
	./holefill

clean:
	rm -f holefill
