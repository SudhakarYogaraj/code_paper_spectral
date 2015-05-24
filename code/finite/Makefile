.SUFFIXES:
.PHONY: all clean problem purge

CXX    = clang++
CXXFLAGS = -O3 -std=c++11 -Wall
target = dev.exec
objects = main.o Problem.o Solver_hmm.o Solver_spectral.o Gaussian_integrator.o toolbox.o io.o tictoc.o

all : $(target)

$(target) : $(objects)
	$(CXX) $(CXXFLAGS) $(objects) -o $(target)

main.o : Problem.hpp Solver_hmm.hpp Solver_spectral.hpp templates.hpp io.hpp structures.hpp
Problem.o : toolbox.hpp
Solver_hmm.o : Problem.hpp toolbox.hpp templates.hpp structures.hpp
Solver_spectral.o : Problem.hpp Gaussian_integrator.hpp toolbox.hpp structures.hpp

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

Problem.cpp : matlab/aux_problem.cpp matlab/build_problem_1d.m
	$(MAKE) -C matlab

problem:
	make -C matlab

clean:
	rm -f $(target) $(objects)

purge:
	rm -f $(target) $(objects)
	make clean -C matlab
