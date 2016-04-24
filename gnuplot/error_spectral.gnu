# Error Plot
reset

# Directory of input and output
input = "./"
out = "./"

# Tics and format
set format y "$10^{%L}$"

# Style for lines and points
# set style line 1 lt 1 lw 1.5 pt 13 ps 1.5 pi 1
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1

unset key

## Degree - error

# Logarithmic scale for y
set logscale y
set logscale x 2

set xlabel "Degree of approximation ($d$)"
set title "Relative error for the homogenized coefficients"

# Data from file
data = "<paste ".input."degree ".input."error"

# EPS-LATEX
set term epslatex
set output out."degree_error.tex"
plot data with linespoints ls 1

# PDF
set term pdf
set output out."_degree_error.pdf"
plot data with linespoints ls 1

## Time - error
set xlabel "Time of computations"
set title "Relative error for the homogenized coefficients"

# Data from file
data = "<paste ".input."time ".input."error"

# EPS-LATEX
set term epslatex
set output out."time_error.tex"
# plot data with points

# PDF
set term pdf
set output out."_time_error.pdf"
# plot data with points
