#!/bin/bash
for (( i=1; i<=100; i++ )); do        
    echo $i
    cp nve.in wat216.prmtop $i
    cd $i
    sander –O –i nve.in –c sample.rst –p wat216.prmtop –r nve.rst –o nve.out 
    cd ..
done
