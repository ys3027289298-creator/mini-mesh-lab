CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
CPPFLAGS += -Iinclude

all: probe test_split

probe: src/probe.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: test_split
	./test_split

test_split: test/test_split.c src/split.c include/mml_split.h
	$(CC) $(CFLAGS) $(CPPFLAGS) test/test_split.c src/split.c -o $@ -lm

diagnose: probe
	./probe

clean:
	rm -f split probe test_split

.PHONY: all test diagnose clean
