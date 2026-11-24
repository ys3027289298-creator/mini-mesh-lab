CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: proxy

proxy: src/proxy.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: proxy
	./proxy

clean:
	rm -f proxy

.PHONY: all test clean
