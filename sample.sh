#!/bin/bash
for (( i=1; i<=100; i++ )); do        
    echo $i
    sander -O -i sample.in -c nvt.rst.4 -p wat216.prmtop -r sample.rst -o sample.out 
    mkdir $i
    cp sample.rst $i/ 
    mv sample.rst nvt.rst.4       
done
