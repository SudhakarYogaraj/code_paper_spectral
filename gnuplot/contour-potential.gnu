reset
f(x,y)=((x - 1*cos(0))**2 + (y - 1*sin(0))**2) * ((x - 1*cos(2 * pi / 3))**2 + (y - 1*sin(2 * pi / 3))**2) * ((x - 1*cos(4 * pi / 3))**2 + (y - 1*sin(4 * pi / 3))**2)
set xrange [-1.4:1.4]
set yrange [-1.4:1.4]
set isosample 250, 250
set table 'test.dat'
splot f(x,y)
unset table

set contour base
set cntrparam level incremental 0, 0.4, 2
unset surface
set table 'cont.dat'
splot f(x,y)
unset table

reset
set xrange [-1.4:1.4]
set yrange [-1.4:1.4]
unset key
set palette rgbformulae 33,13,10

set term pdf
set output "contour-triple-well.pdf"
set term epslatex
set output "contour-triple-well.tex"
set title "Contour lines of the triple-well potential"
set xlabel "$y_0$"
set ylabel "$y_1$"
set cbrange [0:2]
set size square
set xtics -1.2,.4,1.2 out
set ytics -1.2,.4,1.2 out
plot 'test.dat' with image, 'cont.dat' w l lt -1 lw 1.5
