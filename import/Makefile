CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: import

import: src/import.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: import
	./import --selftest

# 生成基准数据并跑三种模式的分段计时
bench: import
	python3 -I tools/gen_obj.py /tmp/bench.obj 250
	./import --legacy /tmp/bench.obj
	./import /tmp/bench.obj
	./import --fast /tmp/bench.obj

clean:
	rm -f import
