CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: collision

collision: src/collision.c src/collision.h
	$(CC) $(CFLAGS) src/collision.c -o $@ -lm

bench: src/bench.c src/collision.c src/collision.h
	$(CC) $(CFLAGS) -DCOLLISION_NO_MAIN src/bench.c src/collision.c -o $@ -lm

bench_ext: src/bench_ext.c src/collision.c src/collision.h
	$(CC) $(CFLAGS) -DCOLLISION_NO_MAIN src/bench_ext.c src/collision.c -o $@ -lm

test: collision
	./collision

clean:
	rm -f collision bench bench_ext
