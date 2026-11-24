CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: attrmap

attrmap: src/attrmap.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: attrmap
	./attrmap

clean:
	rm -f attrmap
