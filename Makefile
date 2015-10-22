# Delete any default suffixes
.SUFFIXES:

# Declare phone targets (always out of date)
.PHONY: all clean clean-dep problem clean-problem

# Compiler and flags
CXX = clang++
CXXFLAGS = -Isrc -O3 -Ofast -ffast-math -std=c++11 -Wall

# Problem file
PRB = src/problem/Problem.cpp

# C++ source files and location of .o files
CPP_FILES := $(wildcard src/*.cpp) $(wildcard src/*/*.cpp) $(PRB)
DEP_FILES := $(subst src,dep, $(CPP_FILES:.cpp=.d))
OBJ_FILES := $(subst src,obj, $(CPP_FILES:.cpp=.o))

# Target executable
TARGET = $(notdir $(shell git rev-parse --abbrev-ref HEAD)).exec

# Program to generate dependencies
MAKEDEPEND = $(CXX) $(CXXFLAGS) -MM -o dep/$*.d $<

# Ensure that directories exist and build executable
all :
	@make --no-print-directory prebuild
	@make --no-print-directory target

prebuild :
	@mkdir -p out $(dir $(DEP_FILES) $(OBJ_FILES))
	cp python/outputs/${ARGS}.cpp $(PRB)

target : $(TARGET)

# Rule fo the target executable
$(TARGET) : $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Include dependencies
-include $(DEP_FILES)

# Create object files from c++ code, and generate dependencies at the same time
obj/%.o : src/%.cpp
	$(MAKEDEPEND)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@sed -i "s#^\([^.]*\.o\)#$(OBJ_DIR)/\1#g" dep/$*.d

problem:
	make -C python

clean:
	rm -f $(TARGET) $(OBJ_FILES) $(DEP_FILES) $(PRB)

clean-problem:
	make clean -C python

clean-all:
	git clean -dxf
