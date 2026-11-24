CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: welded weld_det

welded: src/welded.c
	$(CC) $(CFLAGS) $< -o $@ -lm

weld_det: src/weld_det.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: welded
	./welded

clean:
	rm -f welded weld_det
