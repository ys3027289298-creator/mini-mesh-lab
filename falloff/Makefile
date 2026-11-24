CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: falloff

falloff: src/falloff.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: falloff
	./falloff

clean:
	rm -f falloff
