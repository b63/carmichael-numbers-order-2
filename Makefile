# Program for compiling C++ programs; default ‘g++’.
CXX=g++

# preprocessor flags
override CPPFLAGS:=$(CPPFLAGS) -I ./include

# Extra flags to give to the C++ compiler. 
override CXXFLAGS:=$(CXXFLAGS) -Wall -pedantic -Werror -std=c++14

# Extra flags to give to compilers when they are supposed to invoke the linker, ‘ld’, such as -L. Libraries (-lfoo) should be added to the LDLIBS variable instead. 
override LDFLAGS:=$(LDFLAGS)

# Library flags or names given to compilers when they are supposed to invoke the linker, ‘ld’.
LDLIBS:=-lntl -lgmp -lm -pthread

COMPILE=$(CXX) -c $(CXXFLAGS) $(CPPFLAGS)
LINK=$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) 

%.o : %.cpp
	$(COMPILE) $< -o $@


SOURCES=generate_cprimes.cpp generate_cprimes_order_2.cpp timer.cpp \
		construct_P.cpp counting_factors.cpp
OBJECTS=$(SOURCES:%.cpp=%.o)
BINARIES=generate_cprimes generate_cprimes_order_2 construct_P gen_distributions

.PHONY: all clean

all: $(BINARIES)

$(OBJECTS): include/*

generate_cprimes: generate_cprimes.o timer.o util.o
	$(LINK) $^ -o $@ $(LDLIBS)

generate_cprimes_order_2: generate_cprimes_order_2.o timer.o util.o
	$(LINK) $^ -o $@ $(LDLIBS)

construct_P: construct_P.o timer.o util.o counting_factors.o primality.o
	$(LINK) $^ -o $@ $(LDLIBS)

gen_distributions: gen_distributions.o util.o timer.o
	$(LINK) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(OBJECTS) $(BINARIES)




