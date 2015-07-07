# Delete any default suffixes
.SUFFIXES:

# Declare phone targets (always out of date)
.PHONY: all cpp problem clean clean-problem clean-cpp

# Build the whole project
all :
	make -C python
	make -C cpp

# Build only the cpp part
cpp : 
	make -C cpp

# Build problem (generates Problem.cpp)
problem:
	make -C python

# Clean the whole project
clean:
	make clean -C cpp
	make clean -C python

# Clean cpp files
clean-cpp:
	make clean -C cpp

# Delete problem-specific temporary files
clean-problem:
	make clean -C python

# Clean all non-git files
clean-all:
	git clean -dfx

