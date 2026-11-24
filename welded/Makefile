CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: welded

welded: src/welded.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: welded
	./welded

clean:
	rm -f welded
