CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: skincut

skincut: src/skincut.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: skincut
	./skincut

clean:
	rm -f skincut
