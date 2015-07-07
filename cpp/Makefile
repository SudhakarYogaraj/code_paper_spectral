# Delete any default suffixes
.SUFFIXES:

# Declare phone targets (always out of date)
.PHONY: all clean clean-dep

# Compiler and flags
CXX    = clang++
CXXFLAGS = -O3 -Ofast -fassociative-math -ffast-math -std=c++11 -Wall

# Directories
DEP_DIR = dep
SRC_DIR = src
OBJ_DIR = obj

# C++ source files and location of .o files
CPP_FILES := $(wildcard $(SRC_DIR)/*.cpp)
DEP_FILES := $(addprefix $(DEP_DIR)/,$(notdir $(CPP_FILES:.cpp=.d)))
OBJ_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(CPP_FILES:.cpp=.o)))

# Target executable
TARGET = ../$(notdir $(shell git rev-parse --abbrev-ref HEAD)).exec

# Program to generate dependencies
MAKEDEPEND = $(CXX) $(CXXFLAGS) -MM -o $(DEP_DIR)/$*.d $<

# Ensure that directories exist and build executable
all :
	mkdir -p dep obj
	make target

# Build executable
target : $(TARGET)

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

# Clean the whole project
clean:
	rm -f $(TARGET) $(OBJ_FILES) $(DEP_FILES)

# Delete dependencies
clean-dep:
	rm -f $(DEP_FILES)
