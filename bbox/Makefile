CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Iinclude -I../contract

all: bbox

bbox: src/bbox.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: bbox
	./bbox

clean:
	rm -f bbox
