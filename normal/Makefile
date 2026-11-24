CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Iinclude -I../contract

all: normal

normal: src/normal.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: normal
	./normal

clean:
	rm -f normal
