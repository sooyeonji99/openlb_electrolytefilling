OLB_ROOT := /home/soooo/OpenLB/olb-release
-include $(OLB_ROOT)/config.mk

CXX ?= g++
CXXFLAGS ?= -O3 -std=c++20 -Wall -Wextra

ifeq ($(PARALLEL_MODE),MPI)
PARALLEL_FLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS)
endif

all: nmc_timecost multiphase_smoke

multiphase_smoke: multiphase_smoke.o
	$(CXX) $< -o $@ -lolbcore -L$(OLB_ROOT)/external/lib -lpthread -lz -ltinyxml2 -L$(OLB_ROOT)/build/lib

multiphase_smoke.o: multiphase_smoke.cpp
	$(CXX) $(CXXFLAGS) $(PARALLEL_FLAGS) -DOLB_VERSION=\"NMC_timecost\" -DPLATFORM_CPU_SISD -DDEFAULT_FLOATING_POINT_TYPE=double \
	  -I$(OLB_ROOT)/src -I$(OLB_ROOT)/external/zlib -I$(OLB_ROOT)/external/tinyxml2 -c -o $@ $<

nmc_timecost: main.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f nmc_timecost multiphase_smoke multiphase_smoke.o

.PHONY: all clean
