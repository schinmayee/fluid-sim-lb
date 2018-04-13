ROOT_DIR=..

-include $(ROOT_DIR)/makeinclude/paths.mk
include $(ROOT_DIR)/makeinclude/flags.mk

IFLAGS += -I $(ROOT_DIR)
SHARED_FLAGS = -shared -fPIC

SRC = $(wildcard *.cc)
OBJ = $(patsubst %.cc,build/%.o,$(notdir $(SRC)))

# default build common library
default: common

# shared library
common = build/libcommon.a
common: $(OBJ)
	ar rvs $(common) $^
	#$(CXX) $(CXXFLAGS) $(IFLAGS) $(SHARED_FLAGS) -o $(common) $^ $(LFLAGS)

# object files
build/%.o: %.cc
	@echo "Building $(CURDIR)/$@"
	mkdir -p build
	$(CXX) $(CXXFLAGS) -c $(IFLAGS) $(SHARED_FLAGS) -o $@ $< $(LFLAGS)

# clean
.PHONY: clean
clean:
	@echo "Removing object files from $(CURDIR)/build"
	-rm -f build/*

# purge
.PHONY: purge
purge:
	@echo "Removing object files and build dependencies for $(CURDIR)"
	-ls
	-rm -rf build

# dependencies
build/particles.o: particles.cc particles.h primitives.h metagrid.h
build/primitives.o: primitives.cc primitives.h
build/metagrid.o: metagrid.cc metagrid.h primitives.h
build/scalar_grid.o: scalar_grid.cc scalar_grid.h primitives.h metagrid.h
build/vector_grid.o: vector_grid.cc vector_grid.h primitives.h metagrid.h
build/solver.o: solver.cc solver.h definitions.h metagrid.h vector.h matrix.h scalar_grid.h definitions.h primitives.h
build/solver_driver.o: solver_driver.cc solver_driver.h solver.h definitions.h metagrid.h vector.h matrix.h scalar_grid.h definitions.h primitives.h parallel_utils.h
build/vector.o: vector.cc vector.h scalar_grid.h primitives.h metagrid.h
build/matrix.o: matrix.cc matrix.h primitives.h
build/parallel_utils.o: parallel_utils.cc parallel_utils.h
build/distribution.o: distribution.cc distribution.h metagrid.h primitives.h utils.h
build/flip_simulation.o: flip_simulation.cc flip_simulation.h particles.h primitives.h metagrid.h scalar_grid.h vector_grid.h distribution.h definitions.h
build/flip_app.o: flip_app.cc flip_app.h flip_simulation.h particles.h primitives.h metagrid.h scalar_grid.h vector_grid.h distribution.h definitions.h solver_driver.h solver.h
build/flip_app_sim.o: flip_app_sim.cc flip_app.h flip_simulation.h particles.h primitives.h metagrid.h scalar_grid.h vector_grid.h distribution.h definitions.h solver_driver.h solver.h
