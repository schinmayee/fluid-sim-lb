ROOT=$(abspath ../..)

include $(ROOT)/makeinclude/paths.mk
include $(ROOT)/makeinclude/flags.mk

CXXFLAGS+=-fpic 

IFLAGS+=-I $(ROOT)
LCOMMON:=-Wl,-rpath,$(ROOT)/common/build -L$(ROOT)/common/build -lcommon
LFLAGS+=-lrt -lgflags -lglog
LFLAGS+=-L$(CURDIR)/build

SRC = $(wildcard *.cc)
OBJ = $(patsubst %.cc,build/%.o,$(notdir $(SRC)))
TGT = build/convert

default: TGT

# target
TGT: $(SRC)
	mkdir -p  build
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o $(TGT) $^ $(LCOMMON) $(LFLAGS)

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
build/convert: $(SRC)
