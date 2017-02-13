SDSL_DIR=../sdsl-lite
RLZAP_DIR=../rlzap

# In OS X, getrusage() returns maximum resident set size in bytes.
# In Linux, the value is in kilobytes, so this line should be commented out.
#RUSAGE_FLAGS=-DRUSAGE_IN_BYTES

# Compute run and gap measures for the bitvectors in RelativeFM::reportSize().
# This makes reportSize() significantly slower.
#RUN_FLAGS=-DREPORT_RUNS

# Print some additional information. Status info goes to stderr.
VERBOSE_FLAGS=-DVERBOSE_OUTPUT
#VERBOSE_FLAGS=-DVERBOSE_OUTPUT -DVERBOSE_STATUS_INFO

# Hybrid bitvectors are slower, but they can sometimes be smaller.
#VECTOR_FLAGS=-DUSE_HYBRID_BITVECTORS

# Multithreading with OpenMP. No longer compiles without OpenMP support.
# Currently used for RFM construction.
PARALLEL_FLAGS=-fopenmp -D_GLIBCXX_PARALLEL

OTHER_FLAGS=$(RUSAGE_FLAGS) $(RUN_FLAGS) $(VERBOSE_FLAGS) $(VECTOR_FLAGS) $(PARALLEL_FLAGS)

include $(SDSL_DIR)/Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(OTHER_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -I$(RLZAP_DIR)/include -I$(RLZAP_DIR)/ext_libs/boost/include -I$(RLZAP_DIR)/ext_libs/sais/include
LIBOBJS=relative_fm.o utils.o support.o relative_lcp.o
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
OBJS=$(SOURCES:.cpp=.o)
LIBS=-L$(LIB_DIR) -L$(RLZAP_DIR)/build -lsdsl -ldivsufsort -ldivsufsort64 -lrlz_lib
LIBRARY=librfm.a
PROGRAMS=build_bwt align_bwts build_rlcp query_test cst_traverse cst_compare mutate

all: $(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

$(LIBRARY):$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

build_bwt:build_bwt.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

align_bwts:align_bwts.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

build_rlcp:build_rlcp.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

query_test:query_test.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

cst_traverse:cst_traverse.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

cst_compare:cst_compare.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

mutate:mutate.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

clean:
	rm -f $(PROGRAMS) $(LIBRARY) $(OBJS)
