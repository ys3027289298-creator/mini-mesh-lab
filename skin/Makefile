CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: skin

skin: src/skin.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: skin
	./skin

clean:
	rm -f skin
