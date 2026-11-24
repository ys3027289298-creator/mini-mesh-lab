CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: animsync

animsync: src/animsync.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: animsync
	./animsync

clean:
	rm -f animsync
