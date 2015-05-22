.SUFFIXES:

CXX    = clang++
target = dev.exec
CXXFLAGS = -O3 -std=c++11 -Wall

objects = main.o Problem.o Solver_hmm.o Solver_spectral.o Gaussian_integrator.o toolbox.o io.o tictoc.o
$(target) : $(objects)
	$(CXX) $(CXXFLAGS) $(objects) -o $(target)

main.o : Problem.hpp Solver_hmm.hpp Solver_spectral.hpp templates.hpp io.hpp structures.hpp
Problem.o : toolbox.hpp
Solver_hmm.o : Problem.hpp toolbox.hpp templates.hpp structures.hpp
Solver_spectral.o : Problem.hpp Gaussian_integrator.hpp toolbox.hpp structures.hpp

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(target) $(objects)
