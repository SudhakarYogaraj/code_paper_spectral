cc     = g++
target = main.exec
cflags  = -g -std=c++11 -Wall

default : $(target)

$(target) : main.o Problem.o Solver_hmm.o Solver_spectral.o aux.o
	$(cc) $(cflags) main.o Problem.o Solver_hmm.o Solver_spectral.o aux.o -o $(target)

Problem.o : Problem.cpp
	$(cc) $(cflags) -c Problem.cpp

main.o : main.cpp aux.hpp
	$(cc) $(cflags) -c main.cpp

Solver_hmm.o : Solver_hmm.cpp
	$(cc) $(cflags) -c Solver_hmm.cpp

Solver_spectral.o : Solver_spectral.cpp
	$(cc) $(cflags) -c Solver_spectral.cpp

aux.o : aux.cpp
	$(cc) $(cflags) -c aux.cpp

clean:
	rm $(target) *.o
