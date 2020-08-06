# Program for compiling C++ programs; default ‘g++’.
CXX=g++

# preprocessor flags
override CPPFLAGS:=$(CPPFLAGS) -I ./include

# Extra flags to give to the C++ compiler. 
override CXXFLAGS:=$(CXXFLAGS) -Wall -Werror -pedantic -std=c++14 -Ofast -march=native

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
SOURCES+= calc_density.cpp calc_density_batch.cpp
SOURCES+= $(addprefix strategy_1/, nonrigid.cpp gen_nonrigid.cpp pnonrigid.cpp)
SOURCES+= $(addprefix strategy_2/, nonrigid.cpp gen_nonrigid.cpp, all_possible_nonrigid_pairs.cpp)
SOURCES+= benchmarks/bench_construct_P.cpp 
SOURCES+= benchmarks/bench_prodcache.cpp 
SOURCES+= benchmarks/bench_lambda.cpp 
SOURCES+= $(addprefix tests/, test_subsetprod_mod.cpp test_binomial.cpp )


OBJECTS=$(SOURCES:%.cpp=%.o)

BINARIES = generate_cprimes generate_cprimes_order_2 construct_P gen_distributions
BINARIES+= calc_density calc_density_batch
BINARIES+= $(addprefix strategy_1/, gen_nonrigid pnonrigid)
BINARIES+= $(addprefix strategy_2/, gen_nonrigid all_nonrigid_pairs)

BENCHMARKS = $(addprefix benchmarks/, bench_construct_P bench_prodcache)
TESTS = $(addprefix tests/, test_subsetprod_mod test_binomial)

DIR_GUARD = [ -d $(@D) ] || mkdir -p $(@D)

# STATIC PATTERN MATCHING RULES
.SUFFIXES:
SHELL=/usr/bin/bash

${BUILD_DIR}/%.o: ${SRC_DIR}/%.cpp
	@$(DIR_GUARD)
	$(COMPILE) $< -o $@



# automatically generate header dependencies using -MM preprocessor flag
# http://www.gnu.org/software/make/manual/make.html#Automatic-Prerequisites
# match .cpp files in src/ and in src/benchmarks/
${BUILD_DIR}/deps/%.d: ${SRC_DIR}/%.cpp
	@set -e;\
	$(DIR_GUARD);\
	rm -f $@;\
	$(CXX) -MM $(CPPFLAGS) $< > $@.$$$$;\
	sed 's,\($(basename $(notdir $@))\)\.o[ :]*,${BUILD_DIR}/$*\.o $(basename $@)\.d : ,g' < $@.$$$$ > $@;\
	rm -f $@.$$$$;

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

.PHONY: all benchmarks clean

all: $(addprefix ${BUILD_DIR}/,$(BINARIES))

benchmarks: $(addprefix ${BUILD_DIR}/,$(BENCHMARKS))

tests: $(addprefix ${BUILD_DIR}/,$(TESTS))

clean:
	rm -rf ${BUILD_DIR}/*


${BUILD_DIR}/generate_cprimes: $(addprefix ${BUILD_DIR}/,util.o generate_cprimes.o )
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)

${BUILD_DIR}/generate_cprimes_order_2: $(addprefix ${BUILD_DIR}/,generate_cprimes_order_2.o \
	    timer.o util.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)

${BUILD_DIR}/construct_P: $(addprefix ${BUILD_DIR}/,construct_P_main.o construct_P.o \
	    timer.o util.o counting_factors.o primality.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)

${BUILD_DIR}/gen_distributions: $(addprefix ${BUILD_DIR}/,gen_distributions.o util.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)

${BUILD_DIR}/calc_density: $(addprefix ${BUILD_DIR}/,util.o timer.o calc_density.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)

${BUILD_DIR}/calc_density_batch: $(addprefix ${BUILD_DIR}/,util.o timer.o calc_density_batch.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)


# STRATEGY 1 TARGETS
# =============================

S1_DIR = ${BUILD_DIR}/strategy_1

${S1_DIR}/gen_nonrigid: $(addprefix ${BUILD_DIR}/,util.o) \
		$(addprefix ${S1_DIR}/, nonrigid.o gen_nonrigid.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)

${S1_DIR}/pnonrigid: $(addprefix ${BUILD_DIR}/,util.o timer.o) \
		$(addprefix ${S1_DIR}/, nonrigid.o pnonrigid.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)


# STRATEGY 2 TARGETS
# =============================

S2_DIR = ${BUILD_DIR}/strategy_2

${S2_DIR}/gen_nonrigid: $(addprefix ${BUILD_DIR}/, counting_factors.o util.o) \
		$(addprefix ${S2_DIR}/, nonrigid.o gen_nonrigid.o subset_product.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)

${S2_DIR}/all_nonrigid_pairs: $(addprefix ${BUILD_DIR}/,util.o) \
		$(addprefix ${S2_DIR}/, nonrigid.o all_possible_nonrigid_pairs.o)
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)


# TEST TARGETS
# =============================

TEST_DIR = ${BUILD_DIR}/tests

${TEST_DIR}/test_subsetprod_mod: $(addprefix ${BUILD_DIR}/,util.o) \
		$(addprefix ${S2_DIR}/, subset_product.o ) \
		$(addprefix ${TEST_DIR}/, test_subsetprod_mod.o )
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)

${TEST_DIR}/test_binomial: $(addprefix ${BUILD_DIR}/,util.o) \
		$(addprefix ${S2_DIR}/, subset_product.o ) \
		$(addprefix ${TEST_DIR}/, test_binomial.o )
	$(LINK) $(filter-out %.h,$^) -o $@ $(LDLIBS)


# BENCHMARK TARGETS
# =============================

${BUILD_DIR}/benchmarks/bench_construct_P: $(addprefix ${BUILD_DIR}/,timer.o primality.o \
		util.o counting_factors.o construct_P.o) \
		$(addprefix ${BUILD_DIR}/benchmarks/, bench_construct_P.o)
	$(LINK) $(filter-out %.h,$^) -o $@ -lbenchmark -lbenchmark_main $(LDLIBS) 

${BUILD_DIR}/benchmarks/bench_prodcache: $(addprefix ${BUILD_DIR}/, util.o) \
	    $(addprefix ${BUILD_DIR}/benchmarks/, bench_prodcache.o) \
	    $(addprefix ${BUILD_DIR}/strategy_1/, nonrigid.o)
	$(LINK) $(filter-out %.h,$^) -o $@ -lbenchmark -lbenchmark_main $(LDLIBS) 


include $(foreach src,$(SOURCES),$(BUILD_DIR)/deps/$(basename ${src}).d )
#include $(SOURCES:%.cpp=${BUILD_DIR}/deps/%.d)
#include $(BINARIES:%=${BUILD_DIR}/%.dd)

