CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: holefill

holefill: src/holefill.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: holefill
	./holefill

# run the experiment and deposit PPM renders + metrics JSON in artifacts/
artifacts: holefill
	./holefill ../artifacts/holefill

clean:
	rm -f holefill

.PHONY: all test artifacts clean
