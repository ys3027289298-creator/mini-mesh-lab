MODULES := welded normal bbox

.PHONY: all test modules pipeline clean $(MODULES)

all: test

test: modules pipeline
	@echo "==> all module and pipeline checks passed"

modules:
	@for d in $(MODULES); do \
		$(MAKE) -C $$d test || exit 1; \
	done

pipeline:
	$(MAKE) -C pipeline test

clean:
	@for d in $(MODULES); do \
		$(MAKE) -C $$d clean; \
	done
	$(MAKE) -C pipeline clean
