CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: hull

hull: src/hull.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: hull
	./hull

clean:
	rm -f hull
