CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: inertia

inertia: src/inertia.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: inertia
	./inertia

clean:
	rm -f inertia

.PHONY: all test clean
