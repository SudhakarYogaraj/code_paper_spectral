# Error Plot
reset

# Directory of input and output
input = "out/error_spectral/"
out = "pictures/"

# Terminal
set term pdf

# Style for lines and points
set style line 1 lt 1 lw 1.5 pt 13 ps 1.5 pi 1

# Output file
set output out."time_integration.pdf"

# unset key
set xlabel "Time"
set title "Comparison of the exact and numerical solutions"

# Data from file
data_exact = "<paste ".input."time_exact ".input."sol_exact"
data_spectral = "<paste ".input."time_spectral ".input."sol_spectral"
data_hmm = "<paste ".input."time_hmm ".input."sol_hmm"

plot data_exact with lines, data_spectral with lines, data_hmm with lines
