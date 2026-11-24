CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: normal

normal: src/normal.c src/normal.h
	$(CC) $(CFLAGS) src/normal.c -o $@ -lm

# original end-to-end checks, now driven through the split U4 entry point
test: normal unit_tests bit_identity
	./normal
	./unit_tests
	./bit_identity

build/normal_nom.o: src/normal.c src/normal.h
	@mkdir -p build
	$(CC) $(CFLAGS) -Dmain=normal_demo_main -c src/normal.c -o $@

build/normal_ref.o: tests/normal_ref.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c tests/normal_ref.c -o $@

# branch-level tests for the four split units
unit_tests: tests/unit_tests.c build/normal_nom.o
	$(CC) $(CFLAGS) tests/unit_tests.c build/normal_nom.o -o $@ -lm

# bit-for-bit comparison against the pre-split implementation at git HEAD
# (tests/normal_ref.c is a symbol-renamed copy of the original monolith)
bit_identity: tests/bit_identity.c build/normal_ref.o build/normal_nom.o
	$(CC) $(CFLAGS) tests/bit_identity.c build/normal_ref.o build/normal_nom.o -o $@ -lm

clean:
	rm -rf normal unit_tests bit_identity build
