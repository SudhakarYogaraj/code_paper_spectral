# Delete any default suffixes
.SUFFIXES:

# Declare phone targets (always out of date)
.PHONY: all clean clean-dep problems clean-problems

# ---- COMPILER AND FLAGS ----
CXX = clang++
CXXFLAGS = -Isrc -O3 -Ofast -ffast-math -std=c++11 -Wall
LIBS = -larmadillo

# Problem file
PRB = src/problems/problem_${ARG}.cpp

# All sources, corresponding to different problems
ALL_SOURCES := $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)

# For the particular problem
CPP_FILES := $(filter-out src/problems/problem_%.cpp, $(ALL_SOURCES)) $(PRB)
DEP_FILES := $(subst src,dep, $(CPP_FILES:.cpp=.d))
OBJ_FILES := $(subst src,obj, $(CPP_FILES:.cpp=.o))

# Target executable
TARGET = $(notdir $(shell git rev-parse --abbrev-ref HEAD)).exec

# ---- BUILDING THE FILES FOR TESTS ----
ALL_DIRS := $(wildcard src/*)
LIB_DIRS := $(filter-out src/tests src/problems src/main, $(ALL_DIRS))
LIB_CPP  := $(foreach dir, $(LIB_DIRS), $(wildcard $(dir)/*.cpp))
LIB_OBJ  := $(subst src,obj, $(LIB_CPP:.cpp=.o))
LIB_DEP  := $(subst src,dep, $(LIB_CPP:.cpp=.d))

# Problems and tests
PROBLEMS := $(notdir $(basename $(wildcard src/problems/problem_*.cpp)))
TESTS    := $(notdir $(basename $(wildcard src/tests/*.cpp)))

# Executables for the tests
TARGET_TESTS := $(addprefix tests/, $(foreach p, $(PROBLEMS), $(addprefix $(p)/, $(TESTS))))

# Ensure that directories exist and build executable
all : prebuild $(TARGET) $(TARGET_TESTS)

tests/% : $(LIB_OBJ)
	$(CXX) $(LIB) $(CXXFLAGS) $(LIB_OBJ) \
	$(addsuffix .o, $(addprefix obj/tests/, $(notdir $@))) \
	$(addsuffix .o, $(addprefix obj/problems/, $(shell basename $(dir $@)))) \
	-o $@

# ---- CREATE MISSING DIRECTORIES ----
prebuild :
	@rm -f $(TARGET)
	@mkdir -p out $(dir $(DEP_FILES) $(OBJ_FILES) $(TARGET_TESTS))
	@mkdir -p $(addprefix out/,$(notdir $(basename $(wildcard src/tests/*.cpp))))

# ---- CREATE OBJECT FILES ----
# Include dependencies
-include $(DEP_FILES)

# Program to generate dependencies
MAKEDEPEND = $(CXX) $(CXXFLAGS) -MM -o dep/$*.d $<

# Create object files from c++ code, and generate dependencies at the same time
obj/%.o : src/%.cpp
	$(MAKEDEPEND)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@sed -i "s#^\([^.]*\.o\)#obj/$*.o#g" dep/$*.d

# ---- CREATE MAIN EXECTUABLE ----
# Rule fo the target executable
$(TARGET) : $(OBJ_FILES)
	$(CXX) $(LIBS) $(CXXFLAGS) $^ -o $@

# ---- CREATE TEST EXECUTABLES ----
tests/% : $(OBJ_FILES)

# ---- CREATE PROBLEM FILES ----
problems:
	make -C python

# ---- CLEAN CPP, PROBLEMS, OR ALL GITIGNORE ----
clean:
	rm -rf $(TARGET) obj dep

clean-problems:
	make clean -C python

clean-all:
	git clean -dxf
