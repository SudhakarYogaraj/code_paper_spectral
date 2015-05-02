.SUFFIXES:

cxx    = g++
target = dev.exec
cxxflags = -O3 -std=c++11 -Wall

objects = main.o Problem.o Solver_hmm.o Solver_spectral.o Gaussian_integrator.o aux.o io.o
$(target) : $(objects)
	$(cxx) $(cxxflags) $(objects) -o $(target)

main.o : Problem.hpp Solver_hmm.hpp Solver_spectral.hpp templates.hpp io.hpp
Problem.o : aux.hpp
Solver_hmm.o : Problem.hpp aux.hpp templates.hpp
Solver_spectral.o : Problem.hpp Gaussian_integrator.hpp aux.hpp

%.o: %.cpp %.hpp
	$(cxx) $(cxxflags) -c $< -o $@

.PHONY: clean
clean:
	rm $(target) $(objects)
