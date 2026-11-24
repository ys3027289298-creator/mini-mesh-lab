CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: welded

welded: src/welded.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: welded
	./welded

# ASan/UBSan run of the same self-check (includes the out-of-range index case)
sanitize: src/welded.c
	$(CC) -O1 -g -std=c11 -Wall -Wextra -fsanitize=address,undefined -fno-omit-frame-pointer $< -o welded_san -lm
	ASAN_OPTIONS=detect_leaks=0 ./welded_san

clean:
	rm -f welded welded_san

.PHONY: all test sanitize clean
