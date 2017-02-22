# Time integration plot
reset

# Directory of input and output
input = "./"
out = "./"

set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1
set title "Error in time ($\\left( \\mathbb E \\, \\sup_{0\\, \\leq \\,t\\, \\leq \\,T} \\, |X(t) \\,-\\,X_d(t)|^2 \\right)^{1/2}$)"
set xlabel "Degree of approximation ($d$)"
set format y "$10^{%L}$"
set logscale y
# set logscale x 2
unset key

# Data from file
data = "<paste ".input."degree ".input."error"

# Fit
f(x) = a*x + b
fit f(x) data using 1:(log(sqrt($2))) via a,b

set term pdf
set output out."time_integration.pdf"
plot data using 1:(sqrt($2)) with points ls 1, exp(f(x)) lw 2

set term epslatex
set output out."time_error.tex"
plot data using 1:(sqrt($2)) with points ls 1, exp(f(x)) lw 2
