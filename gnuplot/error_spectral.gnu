# Error Plot
reset

# Directory of input and output
input = "out/error_spectral/"
out = "pictures/"

# Terminal
set term pdf

# Style for lines and points
set style line 1 lt 1 lw 1.5 pt 13 ps 1.5 pi 1

unset key


## Degree - error
set output out."degree_error.pdf"

# Logarithmic scale for y
set logscale y
set logscale x

set xlabel "Degree of approximation"
set title "Relative error for the homogenized coefficients"

# Data from file
data = "<paste ".input."degree ".input."error"

plot data with points


## Time - error
set output out."time_error.pdf"

# unset key
set xlabel "Time of computations"
set title "Relative error for the homogenized coefficients"

# Data from file
data = "<paste ".input."time ".input."error"

plot data with points
