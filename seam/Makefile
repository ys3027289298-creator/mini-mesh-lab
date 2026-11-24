CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

OUT ?= ../artifacts/seam

all: seam_test seam_experiment

seam_test: src/seam_test.c src/seam.c src/seam.h
	$(CC) $(CFLAGS) src/seam_test.c src/seam.c -o seam_test -lm

seam_experiment: src/seam_experiment.c src/seam.c src/seam.h
	$(CC) $(CFLAGS) src/seam_experiment.c src/seam.c -o seam_experiment -lm

test: seam_test
	./seam_test

run: seam_experiment
	mkdir -p $(OUT)
	./seam_experiment $(OUT)

plots:
	python3 plot/make_plots.py $(OUT)

clean:
	rm -f seam_test seam_experiment

.PHONY: all test run plots clean
