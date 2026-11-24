CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: interrupt

interrupt: src/interrupt.c
	$(CC) $(CFLAGS) $< -o $@ -lm

# Fault-injection builds: each breaks the atomicity of exactly one step.
#   fault1 -> weld face compaction (in-place remap, degenerate faces visible)
#   fault2 -> normal finalize (per-vertex publish, half-normalized state)
#   fault3 -> bbox (wrong vertex order, box does not match processed prefix)
fault: fault1 fault2 fault3

fault1: src/interrupt.c
	$(CC) $(CFLAGS) -DFAULT_INJECT=1 $< -o $@ -lm
fault2: src/interrupt.c
	$(CC) $(CFLAGS) -DFAULT_INJECT=2 $< -o $@ -lm
fault3: src/interrupt.c
	$(CC) $(CFLAGS) -DFAULT_INJECT=3 $< -o $@ -lm

test: interrupt
	./interrupt

fault-test: fault
	@for f in fault1 fault2 fault3; do \
	  echo "--- $$f (expect a named failure, exit 0) ---"; ./$$f || exit 1; \
	done

clean:
	rm -f interrupt fault1 fault2 fault3
