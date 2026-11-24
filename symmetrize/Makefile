CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: symmetrize

symmetrize: src/symmetrize.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: symmetrize
	./symmetrize

json: symmetrize
	./symmetrize --json

clean:
	rm -f symmetrize

.PHONY: all test json clean
