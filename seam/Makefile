CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: seam

# textual include of production normal.c; -Dmain renames only its test main,
# seam.c restores the token with #undef main before defining its own main
seam: src/seam.c
	mkdir -p out
	$(CC) $(CFLAGS) -Dmain=normal_demo_main $< -o $@ -lm

test: seam
	./seam

plot: seam
	./seam
	python3 -I tools/plot.py out/lum.csv out/lum.svg

clean:
	rm -f seam
	rm -rf out
