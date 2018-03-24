#!/bin/bash

f="HH-distance.dat"
v="H-H"
list=(`awk '{print $2}' $f | sort -n`)
min=${list[0]}
max=${list[${#list[@]}-1]}

# write the file for gnuplot
echo "set term png enhanced

set output 'dist_${v}.png'
set xlabel '${v} (Angstrom)'
set ylabel 'Frequencies'
min=$min
max=$max
bin_width=(max-min)/100.;     # set bin width    
bin_num(x)=floor(x/bin_width)  
rounded(x)=bin_width*(bin_num(x) + 0.5)
plot '${f}' u (rounded(\$2)):(1) t '' smooth frequency with lines
" |  gnuplot
