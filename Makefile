# Delete any default suffixes
.SUFFIXES:

# Declare phone targets (always out of date)
.PHONY: all clean clean-dep clean-problem clean-all problem

# Compiler and flags
CXX    = clang++
CXXFLAGS = -O3 -std=c++11 -Wall

# Directories
DEP_DIR = dep
SRC_DIR = src
OBJ_DIR = obj

# C++ source files and location of .o files
CPP_FILES := $(wildcard $(SRC_DIR)/*.cpp)
DEP_FILES := $(addprefix $(DEP_DIR)/,$(notdir $(CPP_FILES:.cpp=.d)))
OBJ_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(CPP_FILES:.cpp=.o)))

# Target executable
TARGET = dev.exec

# Program to generate dependencies
MAKEDEPEND = $(CXX) $(CXXFLAGS) -MM -o $(DEP_DIR)/$*.d $<


all : $(TARGET)

# Rule fo the target executable
$(TARGET) : $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Include dependencies
-include $(DEP_FILES)

# Create object files from c++ code, and generate dependencies at the same time
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(MAKEDEPEND)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@sed -i "s#^\([^.]*\.o\)#$(OBJ_DIR)/\1#g" $(DEP_DIR)/$*.d

# Build problem (generates Problem.cpp)
problem:
	make -C matlab

# Detele oject files and executable
clean:
	rm -f $(TARGET) $(OBJ_FILES)

# Delete dependencies
clean-dep: 
	rm -f $(DEP_FILES)

# Delete problem-specific temporary files
clean-problem:
	make clean -C matlab

# Delete all temporary files
clean-all:
	rm -f $(TARGET) $(OBJ_FILES)
	rm -f $(DEP_FILES)
	make clean -C matlab
