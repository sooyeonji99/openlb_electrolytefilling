OLB_ROOT := /home/soooo/OpenLB/olb-release
-include $(OLB_ROOT)/config.mk

CXX ?= g++
CXXFLAGS ?= -O3 -std=c++20 -Wall -Wextra

ifeq ($(PARALLEL_MODE),MPI)
PARALLEL_FLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS)
endif

all: main geometry_validate

main: main.o
	$(CXX) $< -o $@ -lolbcore -L$(OLB_ROOT)/external/lib -lpthread -lz -ltinyxml2 -L$(OLB_ROOT)/build/lib

main.o: main.cpp geometry.h physical_range_coupling.h phase_projection.h phase_outlet.h gravity_excitation.h
	$(CXX) $(CXXFLAGS) $(PARALLEL_FLAGS) -DOLB_VERSION=\"NMC_gravity_scale\" -DPLATFORM_CPU_SISD -DDEFAULT_FLOATING_POINT_TYPE=double \
	  -I$(OLB_ROOT)/src -I$(OLB_ROOT)/external/zlib -I$(OLB_ROOT)/external/tinyxml2 -c -o $@ $<

geometry_validate: geometry_validate.cpp geometry.h
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f geometry_validate main main.o

.PHONY: all clean
