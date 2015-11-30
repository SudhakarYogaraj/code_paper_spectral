# Delete any default suffixes
.SUFFIXES:

# Declare phone targets (always out of date)
.PHONY: all clean clean-dep problems clean-problems

# ---- COMPILER AND FLAGS ----
CXX = clang++
CXXFLAGS = -Isrc -O3 -Ofast -ffast-math -std=c++11 -Wall
LIBS = -larmadillo

# ---- BUILDING LIST OF TEST AND LIB FILES ----

# For the particular problem
ALL_CPP := $(wildcard src/*.cpp) $(wildcard src/*/*.cpp)
ALL_DEP := $(subst src,dep, $(ALL_CPP:.cpp=.d))
ALL_OBJ := $(subst src,obj, $(ALL_CPP:.cpp=.o))

LIB_DIRS := $(filter-out src/tests src/problems src/main, $(wildcard src/*))
LIB_CPP  := $(foreach dir, $(LIB_DIRS), $(wildcard $(dir)/*.cpp))
LIB_OBJ  := $(subst src,obj, $(LIB_CPP:.cpp=.o))

# Problems and tests
PROBLEMS := $(notdir $(basename $(wildcard src/problems/problem_*.cpp)))
TESTS    := $(notdir $(basename $(wildcard src/tests/*.cpp)))
TARGETS  := $(addprefix tests/, $(foreach p, $(PROBLEMS), $(addprefix $(p)/, $(TESTS))))

# Target executable, built from main
TARGET = $(notdir $(shell git rev-parse --abbrev-ref HEAD)).exec

# ---- MAIN BUILD ----
all : prebuild $(TARGET) $(TARGETS)

# ---- PREBUILD: CREATE MISSING DIRECTORIES ----
prebuild :
	@mkdir -p out $(dir $(ALL_DEP) $(ALL_OBJ) $(TARGETS))
	@mkdir -p $(addprefix out/,$(notdir $(basename $(wildcard src/tests/*.cpp))))

# ---- CREATE OBJECT FILES ----
# Include dependencies
-include $(ALL_DEP)

# Program to generate dependencies
MAKEDEPEND = $(CXX) $(CXXFLAGS) -MM -o dep/$*.d $<

# Create object files from c++ code, and generate dependencies at the same time
obj/%.o : src/%.cpp
	$(MAKEDEPEND)
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@sed -i "s#^\([^.]*\.o\)#obj/$*.o#g" dep/$*.d

# ---- BUILDING MAIN ----
$(TARGET) : $(LIB_OBJ) obj/main/main.o obj/problems/problem_${ARG}.o
	$(CXX) $(LIBS) $(CXXFLAGS) $^ -o $@

# ---- BUILDING TESTS ----
tests/% : $(ALL_OBJ)
	$(CXX) $(LIBS) $(CXXFLAGS) $(LIB_OBJ) \
	$(addsuffix .o, $(addprefix obj/tests/, $(notdir $@))) \
	$(addsuffix .o, $(addprefix obj/problems/, $(shell basename $(dir $@)))) \
	-o $@

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
