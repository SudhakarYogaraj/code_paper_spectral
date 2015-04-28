.SUFFIXES:

cxx    = g++
target = multiscale.exec
cxxflags = -O3 -std=c++11 -Wall

objects = main.o Problem.o Solver_hmm.o Solver_spectral.o aux.o
$(target) : $(objects)
	$(cxx) $(cxxflags) $(objects) -o $(target)

main.o : Solver_hmm.hpp Solver_spectral.hpp
Problem.o : aux.hpp 
Solver_hmm.o : Problem.hpp aux.hpp
Solver_spectral.o : Problem.hpp aux.hpp

%.o: %.cpp %.hpp
	$(cxx) $(cxxflags) -c $< -o $@

.PHONY: clean
clean:
	rm $(target) $(objects)
