CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Iinclude -I../contract

all: welded

welded: src/welded.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: welded
	./welded

clean:
	rm -f welded
