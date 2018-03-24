#!/bin/bash
sander -O -i nvt.in -p wat216.prmtop -c npt6.rst  -o nvt.out.1 -r nvt.rst.1 -x mdcrd7
sander -O -i nvt.in -p wat216.prmtop -c nvt.rst.1 -o nvt.out.2 -r nvt.rst.2 -x mdcrd8
sander -O -i nvt.in -p wat216.prmtop -c nvt.rst.2 -o nvt.out.3 -r nvt.rst.3 -x mdcrd9
sander -O -i nvt.in -p wat216.prmtop -c nvt.rst.3 -o nvt.out.4 -r nvt.rst.4 -x mdcrd10
