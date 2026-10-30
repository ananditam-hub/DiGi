# Build DiGi's C++ tools into bin/
#   make         node k-mer and KMC k-mer set samplers (no Gurobi needed)
#   make ilp     ILP solver (needs GUROBI_HOME, e.g. .../gurobi1300/linux64)
#   make clean   remove bin/

CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17

# Installed Gurobi library version, e.g. libgurobi130.so -> gurobi130
GRB_LIB = $(patsubst lib%.so,%,$(notdir $(firstword $(wildcard $(GUROBI_HOME)/lib/libgurobi[0-9]*.so))))

all: bin/node_kmers_gen bin/kmer_set_sampling

ilp: bin/dip_ilp

bin/node_kmers_gen: src/graph/node_kmers_gen.cpp | bin
	$(CXX) $(CXXFLAGS) $< -o $@

bin/kmer_set_sampling: src/kmer_sampling/kmer_set_sampling.cpp | bin
	$(CXX) $(CXXFLAGS) $< -o $@

bin/dip_ilp: src/ilp/dip_main_recomb.cpp | bin
	@test -n "$(GUROBI_HOME)" || { echo "set GUROBI_HOME to your Gurobi install" >&2; exit 1; }
	@test -n "$(GRB_LIB)" || { echo "no libgurobi*.so in $(GUROBI_HOME)/lib" >&2; exit 1; }
	$(CXX) $(CXXFLAGS) $< -o $@ -I$(GUROBI_HOME)/include -L$(GUROBI_HOME)/lib \
	  -Wl,-rpath,$(GUROBI_HOME)/lib -lgurobi_c++ -l$(GRB_LIB) -lstdc++fs

bin:
	mkdir -p bin

clean:
	rm -rf bin

.PHONY: all ilp clean
