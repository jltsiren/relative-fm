SDSL_DIR=../sdsl-lite

# In OS X, getrusage() returns maximum resident set size in bytes.
# In Linux, the value is in kilobytes, so this line should be commented out.
#RUSAGE_FLAGS=-DRUSAGE_IN_BYTES

# Return information on the number of runs of 1-bits in RelativeFM::reportSize().
#RUN_FLAGS=-DREPORT_RUNS

# Print some additional information.
VERBOSE_FLAGS=-DVERBOSE_OUTPUT

# Sparse bitvectors are smaller but slower.
#VECTOR_FLAGS=-DUSE_SPARSE_BITVECTORS

# RRR bitvectors make the WT smaller and slower.
#WT_FLAGS=-DUSE_RRR_WT

OTHER_FLAGS=$(RUSAGE_FLAGS) $(RUN_FLAGS) $(VERBOSE_FLAGS) $(VECTOR_FLAGS) $(WT_FLAGS)

include $(SDSL_DIR)/Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(OTHER_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR)
RFMOBJS=relative_fm.o utils.o
SOURCES=$(wildcard *.cpp)
HEADERS=$(wildcard *.h)
OBJS=$(SOURCES:.cpp=.o)
LIBS=-L$(LIB_DIR) $(RFMOBJS) -lsdsl -ldivsufsort -ldivsufsort64
PROGRAMS=align_bwts build_bwt query_test

all: $(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

align_bwts:align_bwts.o $(RFMOBJS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

build_bwt:build_bwt.o $(RFMOBJS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

query_test:query_test.o $(RFMOBJS)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)

package:
	mkdir relative-fm
	cp $(SOURCES) $(HEADERS) Makefile targz relative-fm.tex relative-fm
	./targz relative-fm
	rm relative-fm/*
	rmdir relative-fm

clean:
	rm -f $(PROGRAMS) $(OBJS)
	rm -f relative-fm.aux relative-fm.log relative-fm.pdf
