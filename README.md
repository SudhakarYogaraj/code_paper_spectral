This repository contains the code used for the numerical experiments presented in our paper entitled *Spectral methods for multiscale stochastic differential equations*.

Our implementation uses the C++ linear algebra library [Armadillo](http://arma.sourceforge.net/), which is packaged for most of the popular linux distributions. To build the code, follow these steps:

- First generate the C++ problem files from their python version. For this, enter the directory */python* and execute `make all`. This will create the C++ problem files in the directory */src/problems*.

- In the root folder, run `make tests` to generate executables that will allow running the tests on the different problems. The tests executables are stored in */tests/problem_name/test_name*.

- Run the tests. In the paper, the tests we use are *time_integration*, to two trajectories, and *error_spectral*, to compare the coefficients at a given point of the slow process.

- Make the figures using the scripts in */gnuplot*.

For question, please feel free to contact me.
