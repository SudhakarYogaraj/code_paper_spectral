# Delete any default suffixes
.SUFFIXES:

# Declare phone targets (always out of date)
.PHONY: all clean clean-dep problem clean-problem

# Compiler and flags
CXX    = clang++
CXXFLAGS = -O3 -Ofast -fassociative-math -ffast-math -std=c++11 -Wall

# C++ source files and location of .o files
CPP_FILES := $(wildcard src/*.cpp)
DEP_FILES := $(addprefix dep/,$(notdir $(CPP_FILES:.cpp=.d)))
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# Target executable
TARGET = $(notdir $(shell git rev-parse --abbrev-ref HEAD)).exec

# Program to generate dependencies
MAKEDEPEND = $(CXX) $(CXXFLAGS) -MM -o dep/$*.d $<

# Ensure that directories exist and build executable
all : prebuild $(TARGET)

prebuild :
	mkdir -p dep obj out
	cp python/outputs/${ARGS}.cpp src/Problem.cpp

# Rule fo the target executable
$(TARGET) : $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Include dependencies
-include $(DEP_FILES)

# Create object files from c++ code, and generate dependencies at the same time
obj/%.o: src/%.cpp
	$(MAKEDEPEND)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@sed -i "s#^\([^.]*\.o\)#$(OBJ_DIR)/\1#g" dep/$*.d

clean:
	rm -f $(TARGET) $(OBJ_FILES) $(DEP_FILES)

## Handling of user input ##

problem:
	make -C python

clean-problem:
	make clean -C python
