# compiler
CXX = g++

CXXFLAGS += -std=c++14

OPTIMIZE := -Ofast -DNDEBUG

# check what is available
has_blosc := no
ifneq (,$(and $(BLOSC_LIB_DIR),$(BLOSC_INCL_DIR),$(BLOSC_LIB)))
    has_blosc := yes
endif
has_log4cplus := no
ifneq (,$(and $(LOG4CPLUS_LIB_DIR),$(LOG4CPLUS_INCL_DIR),$(LOG4CPLUS_LIB)))
    has_log4cplus := yes
endif

# compiler flags
CXXFLAGS += -pthread $(OPTIMIZE)
ifeq (yes,$(has_blosc))
    CXXFLAGS += -DOPENVDB_USE_BLOSC
endif
ifeq (yes,$(has_log4cplus))
    CXXFLAGS += -DOPENVDB_USE_LOG4CPLUS
endif
ifeq (2,$(strip $(abi)))
    CXXFLAGS += -DOPENVDB_2_ABI_COMPATIBLE
endif
ifneq (2,$(strip $(GLFW_MAJOR_VERSION)))
    CXXFLAGS += -DOPENVDB_USE_GLFW_3
endif

# include flags
IFLAGS := -I $(OPENVDB_INCL_DIR) -isystem $(BOOST_INCL_DIR) \
  -isystem $(ILMBASE_INCL_DIR) -isystem $(TBB_INCL_DIR)
ifeq (yes,$(has_blosc))
    IFLAGS += -isystem $(BLOSC_INCL_DIR)
endif
ifeq (yes,$(has_log4cplus))
    IFLAGS += -isystem $(LOG4CPLUS_INCL_DIR)
endif
IFLAGS += -isystem $(CEREAL_INCL_DIR)

IFLAGS+=-I $(CANARY_INCL_DIR)
IFLAGS+=-I $(CANARY_DEP_INCL_DIR)

# lib flags
LIBS_RPATH := \
  -Wl,-rpath,$(OPENVDB_LIB_DIR) -L$(OPENVDB_LIB_DIR) $(OPENVDB_LIB)\
  -ldl -lm -lz \
  -Wl,-rpath,$(ILMBASE_LIB_DIR) -L$(ILMBASE_LIB_DIR) $(HALF_LIB) \
  -Wl,-rpath,$(TBB_LIB_DIR) -L$(TBB_LIB_DIR) $(TBB_LIB) \
  -Wl,-rpath,$(BOOST_LIB_DIR) -L$(BOOST_LIB_DIR) $(BOOST_LIB)
ifeq (yes,$(has_blosc))
    LIBS += -L$(BLOSC_LIB_DIR) $(BLOSC_LIB)
    LIBS_RPATH += -Wl,-rpath,$(BLOSC_LIB_DIR) -L$(BLOSC_LIB_DIR) $(BLOSC_LIB)
endif
ifeq (yes,$(has_log4cplus))
    LIBS += -L$(LOG4CPLUS_LIB_DIR) $(LOG4CPLUS_LIB)
    LIBS_RPATH += -Wl,-rpath,$(LOG4CPLUS_LIB_DIR) -L$(LOG4CPLUS_LIB_DIR) $(LOG4CPLUS_LIB)
endif
ifneq (,$(strip $(CONCURRENT_MALLOC_LIB)))
ifneq (,$(strip $(CONCURRENT_MALLOC_LIB_DIR)))
    LIBS_RPATH += -Wl,-rpath,$(CONCURRENT_MALLOC_LIB_DIR) -L$(CONCURRENT_MALLOC_LIB_DIR)
endif
endif
ifdef LINUX
    LIBS += -lrt
    LIBS_RPATH += -lrt
endif
LFLAGS += $(LIBS_RPATH) -L$(CURDIR)
