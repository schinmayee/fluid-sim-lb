# OpenVDB
OPENVDB_INCL_DIR := $(HOME)/installations/openvdb_v4_0/include
OPENVDB_LIB_DIR := $(HOME)/installations/openvdb_v4_0/lib
OPENVDB_LIB := -lopenvdb

# The parent directory of the boost/ header directory
BOOST_INCL_DIR := $(HOME)/installations/boost
# The directory containing libboost_iostreams, libboost_system, etc.
BOOST_LIB_DIR := $(HOME)/installations/boost/stage/lib
BOOST_LIB := -lboost_iostreams
BOOST_LIB += -lboost_program_options -lboost_filesystem
BOOST_LIB += -lboost_system

# The parent directory of the OpenEXR/ header directory
EXR_INCL_DIR := /usr/local/include
# The directory containing IlmImf
EXR_LIB_DIR := /usr/local/lib
EXR_LIB := -lIlmImf

# The parent directory of the OpenEXR/ header directory (which contains half.h)
ILMBASE_INCL_DIR := $(EXR_INCL_DIR)
# The directory containing libIlmThread, libIlmThread, libHalf etc.
ILMBASE_LIB_DIR := $(EXR_LIB_DIR)
ILMBASE_LIB := -lIlmThread -lIex -lImath
HALF_LIB := -lHalf

# The parent directory of the tbb/ header directory
TBB_INCL_DIR := /usr/local/include/tbb
# The directory containing libtbb
TBB_LIB_DIR := /usr/lib/x86_64-linux-gnu
TBB_LIB := -ltbb

# The parent directory of the blosc.h header
# (leave blank if Blosc is unavailable)
BLOSC_INCL_DIR :=
# The directory containing libblosc
BLOSC_LIB_DIR :=
BLOSC_LIB :=

# A scalable, concurrent malloc replacement library
# such as jemalloc (included in the Houdini HDK) or TBB malloc
# (leave blank if unavailable)
CONCURRENT_MALLOC_LIB := -ltbbmalloc_proxy -ltbbmalloc
# The directory containing the malloc replacement library
CONCURRENT_MALLOC_LIB_DIR := /usr/local/intel/composer_xe_2013.5.192/tbb/lib/mic

# The parent directory of the cppunit/ header directory
# (leave blank if CppUnit is unavailable)
CPPUNIT_INCL_DIR := /usr/include
# The directory containing libcppunit
CPPUNIT_LIB_DIR := /usr/lib/x86_64-linux-gnu/
CPPUNIT_LIB := -lcppunit

# The parent directory of the log4cplus/ header directory
# (leave blank if log4cplus is unavailable)
LOG4CPLUS_INCL_DIR := /usr/include
# The directory containing liblog4cplus
LOG4CPLUS_LIB_DIR := /usr/lib/x86_64-linux-gnu/
LOG4CPLUS_LIB := -llog4cplus

# The directory containing glfw.h
# (leave blank if GLFW is unavailable)
GLFW_INCL_DIR := /usr/local/include/
# The directory containing libglfw
GLFW_LIB_DIR := /usr/local/lib
GLFW_LIB := -lglfw3 -lGL -lm -lXrandr -lXi -lX11 -lXxf86vm -lpthread
# The major version number of the GLFW library
# (header filenames changed between GLFW 2 and 3, so this must be specified explicitly)
GLFW_MAJOR_VERSION := 3

# Tetlib
TET_LIB_INCL_DIR := $(HOME)/amr/libs/tetgen
TET_LIB_DIR := $(HOME)/amr/libs/tetgen
TET_LIB := -ltet

# Cereal
CEREAL_INCL_DIR := $(HOME)/research/runtimes/canary/build/dependency/include

# Canary
CANARY := $(HOME)/research/runtimes/canary
CANARY_INCL_DIR := $(CANARY)/include 
CANARY_DEP_INCL_DIR := $(CANARY)/build/dependency/include
