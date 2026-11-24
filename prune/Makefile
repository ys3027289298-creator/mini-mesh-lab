CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: run

run: main.c mesh.c mesh.h
	$(CC) $(CFLAGS) main.c mesh.c -o run -lm

clean:
	rm -f run

.PHONY: all run clean
