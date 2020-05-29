# Program for compiling C++ programs; default ‘g++’.
CXX=g++

# preprocessor flags
CPPFLAGS=

# Extra flags to give to the C++ compiler. 
CXXFLAGS=-Wall -pedantic -Werror -std=c++11

# Extra flags to give to compilers when they are supposed to invoke the linker, ‘ld’, such as -L. Libraries (-lfoo) should be added to the LDLIBS variable instead. 
LDFLAGS=

# Library flags or names given to compilers when they are supposed to invoke the linker, ‘ld’.
LDLIBS=-lntl -lgmp -lm -pthread

COMPILE=$(CXX) -c $(CXXFLAGS) $(CPPFLAGS)
LINK=$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) 

%.o : %.cpp
	$(COMPILE) $< -o $@


SOURCES=generate_cprimes.cpp generate_cprimes_order_2.cpp timer.cpp
OBJECTS=$(SOURCES:%.cpp=%.o)
BINARIES=generate_cprimes generate_cprimes_order_2

.PHONY: all clean

all: $(BINARIES)

$(OBJECTS): timer.h

generate_cprimes: generate_cprimes.o timer.o
	$(LINK) $^ -o $@ $(LDLIBS)

generate_cprimes_order_2: generate_cprimes_order_2.o timer.o
	$(LINK) $^ -o $@ $(LDLIBS)

clean:
	rm $(OBJECTS) $(BINARIES)




