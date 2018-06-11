ROOT_DIR=..

-include $(ROOT_DIR)/makeinclude/paths.mk
include $(ROOT_DIR)/makeinclude/flags.mk

IFLAGS += -I $(ROOT_DIR)
SHARED_FLAGS = -shared -fPIC

HEADERS = $(wildcard *.h)
HEADERS += $(wildcard $(CANARY_INCL_DIR)/*.h)
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
build/%.o: %.cc $(HEADERS)
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
