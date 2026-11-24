CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
CPPFLAGS += -I../contract \
            -I../welded/include \
            -I../normal/include \
            -I../bbox/include
LDLIBS ?= -lm

MODULE_SRCS := ../welded/src/welded.c \
               ../normal/src/normal.c \
               ../bbox/src/bbox.c
CONTRACT_HEADERS := ../contract/mml_contract.h \
                    $(wildcard ../welded/include/*.h) \
                    $(wildcard ../normal/include/*.h) \
                    $(wildcard ../bbox/include/*.h)

.PHONY: all test clean

all: minimal_chain evolution_check

minimal_chain: minimal_chain.c $(MODULE_SRCS) $(CONTRACT_HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -DMML_NO_SELFTEST \
		minimal_chain.c $(MODULE_SRCS) -o $@ $(LDLIBS)

evolution_check: evolution_check.c $(MODULE_SRCS) $(CONTRACT_HEADERS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -DMML_NO_SELFTEST \
		evolution_check.c $(MODULE_SRCS) -o $@ $(LDLIBS)

test: minimal_chain evolution_check
	./minimal_chain
	./evolution_check

clean:
	rm -f minimal_chain evolution_check
