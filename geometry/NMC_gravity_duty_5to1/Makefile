OLB_ROOT ?= $(abspath ../../..)
-include $(OLB_ROOT)/config.mk

# run_mpi.sh requires an MPI build, even when the root config is serial.
CXX := mpic++
PARALLEL_MODE := MPI
CXXFLAGS ?= -O3 -std=c++20 -Wall -Wextra

ifeq ($(PARALLEL_MODE),MPI)
PARALLEL_FLAGS := -DPARALLEL_MODE_MPI $(MPIFLAGS)
endif

all: main geometry_validate

# Keep this model's MPI core separate from the root's serial library.
CORE_SOURCES := communication/mpiManager.cpp communication/ompManager.cpp core/olbInit.cpp core/expr.cpp io/ostreamManager.cpp
CORE_OBJECTS := $(addprefix .build/mpi/,$(CORE_SOURCES:.cpp=.o))

.build/mpi/%.o: $(OLB_ROOT)/src/%.cpp Makefile
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(PARALLEL_FLAGS) -DPLATFORM_CPU_SISD -DDEFAULT_FLOATING_POINT_TYPE=double -I$(OLB_ROOT)/src -MMD -MP -c -o $@ $<

.build/mpi/libolbcore.a: $(CORE_OBJECTS)
	$(AR) rcs $@ $^

main: main.o .build/mpi/libolbcore.a
	$(CXX) main.o .build/mpi/libolbcore.a -o $@ -L$(OLB_ROOT)/external/lib -lpthread -lz -ltinyxml2

main.o: main.cpp geometry.h physical_range_coupling.h phase_projection.h phase_outlet.h gravity_excitation.h Makefile
	$(CXX) $(CXXFLAGS) $(PARALLEL_FLAGS) -DOLB_VERSION=\"NMC_gravity_duty_5to1\" -DPLATFORM_CPU_SISD -DDEFAULT_FLOATING_POINT_TYPE=double \
	  -I$(OLB_ROOT)/src -I$(OLB_ROOT)/external/zlib -I$(OLB_ROOT)/external/tinyxml2 -c -o $@ $<

geometry_validate: geometry_validate.cpp geometry.h
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f geometry_validate main main.o $(CORE_OBJECTS) $(CORE_OBJECTS:.o=.d) .build/mpi/libolbcore.a

-include $(CORE_OBJECTS:.o=.d)

.PHONY: all clean
