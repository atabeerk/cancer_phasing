CXX ?= g++
PYTHON ?= python3
HTSLIB_PREFIX ?= $(CONDA_PREFIX)
UNAME_S := $(shell uname -s)

SRC = main.cpp chromosome_pipeline.cpp pipeline_helpers.cpp snv_analysis.cpp
HEADERS = pipeline.hpp snv_analysis.hpp
BIN = main

CXXFLAGS ?= -std=c++17 -g
CPPFLAGS ?= -D__STDC_LIMIT_MACROS -I$(HTSLIB_PREFIX)/include
LDFLAGS ?= -L$(HTSLIB_PREFIX)/lib
LDLIBS ?= -lhts -lz -lpthread
ifneq ($(UNAME_S),Darwin)
LDLIBS += -lstdc++fs
endif

.PHONY: help build run test test-smoke clean

help:
	@echo "make build       # compile evo/main"
	@echo "make run         # alias for build"
	@echo "make test        # run the bundled BAM/VCF example and check its output"
	@echo "make test-smoke  # lightweight runtime checks"
	@echo "make clean       # remove build artifact"

build: $(BIN)

$(BIN): $(SRC) $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(SRC) $(LDLIBS) -o $(BIN)

run: build

test: build
	$(PYTHON) ../tests/run_toy_test.py

test-smoke: build
	@./$(BIN) --help >/dev/null
	@$(PYTHON) run_genome_by_chr.py --help >/dev/null
	@echo "Smoke test passed."


clean:
	rm -f $(BIN)
