ROOT=$(abspath ../..)

include $(ROOT)/makeinclude/paths.mk
include $(ROOT)/makeinclude/flags.mk

CXXFLAGS+=-fpic 

IFLAGS+=-I $(ROOT)
LCOMMON:=-Wl,-rpath,$(ROOT)/common/build -L$(ROOT)/common/build -lcommon
LFLAGS+=-lrt
LFLAGS+=-L$(CURDIR)/build

COMMON_FILES = $(wildcard $(ROOT)/common/*.h)

SRC = $(wildcard *.cc)
OBJ = $(patsubst %.cc,build/%.o,$(notdir $(SRC)))
TGT = build/water_app.so

default: TGT

# object files
build/%.o: %.cc
	@echo "Building $(CURDIR)/$@"
	mkdir -p build
	$(CXX) $(CXXFLAGS) -c $(IFLAGS) -o $@ $< $(LCOMMON) $(LFLAGS)

# shared library
TGT: $(OBJ)
	mkdir -p  build
	$(CXX) $(CXXFLAGS) $(IFLAGS) -shared -o $(TGT) $^ $(LCOMMON) $(LFLAGS)

# clean
.PHONY: clean
clean:
	@echo "Removing object files from $(CURDIR)/build"
	-rm -f build/*

# dependencies
build/water_app.o: water_app.cc water_app.h $(COMMON_FILES)
