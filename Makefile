UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
    CXX ?= /opt/homebrew/opt/gcc@13/bin/g++-13
else
    CXX ?= g++
endif

CXXFLAGS ?= -std=c++20

CPPFLAGS := $(CPPFLAGS) -I./include -I./include/Utils -I./include/Implementation -I./include/Dynamic -I./include/Compressed
SRCS = $(wildcard ./src/*.cpp)

DOXYGEN ?= doxygen
DOXYFILE ?= Doxyfile

EXEC ?= executable
DOCS ?= documentation
RM ?= rm 

.PHONY = all debug test clean

all:
	@$(CXX) $(CXXFLAGS) -o3 $(CPPFLAGS) $(SRCS) -o $(EXEC)

debug:
#	@$(CXX) $(CXXFLAGS) -fconcepts-diagnostics-depth=3 $(CPPFLAGS) -DDEBUG $(SRCS) -o $(EXEC)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DDEBUG $(SRCS) -o $(EXEC)

test:
#	@$(CXX) $(CXXFLAGS) -fconcepts-diagnostics-depth=3 $(CPPFLAGS) -I./tests -DTEST $(SRCS) tests/tests.cpp -o $(EXEC)
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I./tests/ -DTEST -DDEBUG $(SRCS) tests/tests.cpp -o $(EXEC)

clean:
	@$(RM) *.o *.a
	@$(RM) -f $(EXEC) 
	@$(RM) -rf $(DOCS)

doc:
	@$(DOXYGEN) $(DOXYFILE)

