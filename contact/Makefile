CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: contact

contact: src/contact.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: contact
	./contact

clean:
	rm -f contact
