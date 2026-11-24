CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: tangent

tangent: src/tangent.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: tangent
	./tangent

clean:
	rm -f tangent
