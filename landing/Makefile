CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: landing

landing: src/landing.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: landing
	./landing

clean:
	rm -f landing
