# Program for compiling C++ programs; default ‘g++’.
CXX=g++

# preprocessor flags
override CPPFLAGS:=$(CPPFLAGS) -I ./include

# Extra flags to give to the C++ compiler. 
override CXXFLAGS:=$(CXXFLAGS) -Wall -Werror -pedantic -std=c++14 -O3

# Extra flags to give to compilers when they are supposed to invoke the linker, ‘ld’, such as -L. Libraries (-lfoo) should be added to the LDLIBS variable instead. 
override LDFLAGS:=$(LDFLAGS)

# Library flags or names given to compilers when they are supposed to invoke the linker, ‘ld’.
LDLIBS:=-lpthread -lm -lgmp -lntl

COMPILE=$(CXX) -c $(CXXFLAGS) $(CPPFLAGS)
LINK=$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) 

BUILD_DIR=build
SRC_DIR=src
BENCHMARK_SRC_DIR=${SRC_DIR}/benchmarks

SOURCES = timer.cpp util.cpp primality.cpp  counting_factors.cpp util.cpp
SOURCES+= gen_distributions.cpp
SOURCES+= generate_cprimes.cpp generate_cprimes_order_2.cpp
SOURCES+= construct_P.cpp construct_P_main.cpp
SOURCES+= nonrigid.cpp gen_nonrigid.cpp
SOURCES+= benchmarks/bench_construct_P.cpp 
SOURCES+= benchmarks/bench_prodcache.cpp 
SOURCES+= benchmarks/bench_lambda.cpp 

OBJECTS=$(SOURCES:%.cpp=%.o)
BINARIES=generate_cprimes generate_cprimes_order_2 construct_P gen_distributions \
         gen_nonrigid
BENCHMARKS=bench_construct_P bench_prodcache


DIR_GUARD = [ -d $(@D) ] || mkdir -p $(@D)

# STATIC PATTERN MATCHING RULES
.SUFFIXES:
.ONESHELL:
SHELL=/usr/bin/bash

${BUILD_DIR}/%.o: ${BENCHMARK_SRC_DIR}/%.cpp
	$(COMPILE) $< -o $@

${BUILD_DIR}/%.o: ${SRC_DIR}/%.cpp
	$(COMPILE) $< -o $@



# automatically generate header dependencies using -MM preprocessor flag
# http://www.gnu.org/software/make/manual/make.html#Automatic-Prerequisites
# match .cpp files in src/ and in src/benchmarks/
${BUILD_DIR}/deps/%.d: ${SRC_DIR}/%.cpp
	@set -e
	$(DIR_GUARD)
	rm -f $@
	$(CXX) -MM $(CPPFLAGS) $< > $@.$$$$
	sed 's,\($(basename $@)\)\.o[ :]*,\1.o $(basename $@).d : ,g' < $@.$$$$ > $@
	rm -f $@.$$$$

# do it agains for benchmark source files b/c they are not in src/
# but in src/benchmarks/
# TODO: better way to do this?
${BUILD_DIR}/deps/%.d: ${BENCHMARK_SRC_DIR}/%.cpp
	@set -e
	$(DIR_GUARD)
	rm -f $@
	$(CXX) -MM $(CPPFLAGS) $< > $@.$$$$
	sed 's,\($(basename $@)\)\.o[ :]*,\1.o $(basename $@).d : ,g' < $@.$$$$ > $@
	rm -f $@.$$$$

# automatically generate some of the .o dependencies for binary xx
# * for every included .h file in xx.h, checks if correspoding .o files
#   is in OBJECTS. If it is, then it is added as a prerequisite
${BUILD_DIR}/%.dd: ${BUILD_DIR}/%.d
	@set -e
	$(DIR_GUARD)
	rm -f $@
	printf "$(basename $@) $(basename $@).dd : " > $@
	sed -rE 's,(.*:), ,g; s,\s+,\n,g; s,\\,,g' $(basename $@).d | sed -rEn 's,(\S*/|)(\S+).h,\2.o,p' > $@.$$$$
	while read dep; do
		case "${OBJECTS}" in
		    $$dep\ *|*\ $$dep|*\ $$dep\ *|$$dep)
			printf "${BUILD_DIR}/$$dep " >> $@
			;;
		esac
	done < $@.$$$$
	rm -f $@.$$$$

.PHONY: all benchmarks clean dir

all: $(addprefix ${BUILD_DIR}/,$(BINARIES))

benchmarks: $(addprefix ${BUILD_DIR}/,$(BENCHMARKS))

dir:
	@${DIR_GUARD}
clean:
	rm -rf ${BUILD_DIR}/*

${BUILD_DIR}/generate_cprimes: $(addprefix ${BUILD_DIR}/,util.o generate_cprimes.o )
	$(LINK) $^ -o $@ $(LDLIBS)

${BUILD_DIR}/generate_cprimes_order_2: $(addprefix ${BUILD_DIR}/,generate_cprimes_order_2.o \
            timer.o util.o)
	$(LINK) $^ -o $@ $(LDLIBS)

${BUILD_DIR}/construct_P: $(addprefix ${BUILD_DIR}/,construct_P_main.o construct_P.o \
            timer.o util.o counting_factors.o primality.o)
	$(LINK) $^ -o $@ $(LDLIBS)

${BUILD_DIR}/gen_distributions: $(addprefix ${BUILD_DIR}/,gen_distributions.o util.o)
	$(LINK) $^ -o $@ $(LDLIBS)

${BUILD_DIR}/gen_nonrigid: $(addprefix ${BUILD_DIR}/,util.o nonrigid.o gen_nonrigid.o)
	$(LINK) $^ -o $@ $(LDLIBS)

${BUILD_DIR}/bench_construct_P: $(addprefix ${BUILD_DIR}/,timer.o primality.o util.o counting_factors.o \
	    construct_P.o bench_construct_P.o)
	$(LINK) $^ -o $@ -lbenchmark -lbenchmark_main $(LDLIBS) 

${BUILD_DIR}/bench_prodcache: $(addprefix ${BUILD_DIR}/,bench_prodcache.o nonrigid.o util.o)
	$(LINK) $^ -o $@ -lbenchmark -lbenchmark_main $(LDLIBS) 

include $(foreach src,$(SOURCES),$(BUILD_DIR)/deps/$(basename $(notdir $(src))).d )
#include $(SOURCES:%.cpp=${BUILD_DIR}/deps/%.d)
#include $(BINARIES:%=${BUILD_DIR}/%.dd)

