EXAMPLE = main

OLB_ROOT := /home/soooo/OpenLB/olb-release

SRC = main.cpp \
      geometry.cpp \
      analysis.cpp \
      vtk_writer.cpp

include $(OLB_ROOT)/default.mk

# Generate each dependency file from its corresponding source. The shared
# default rule uses the first source for every .d file in multi-source projects.
%.d: %.cpp
	@$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) $< -MM -MT $(@:.d=.o) >$@

# Track included headers during compilation as well.
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_FLAGS) -MMD -MP -c -o $@ $<
