set term png enhanced
set output 'HH-distance.png'
set xlabel 'nstep'
set ylabel 'H-H (Angstrom)'
plot 'HH-distance.dat' using 1:2 with points title 'HH-distance'
