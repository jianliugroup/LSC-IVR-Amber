#!/bin/bash
# run 1 step using min.rst as coordinate file in the groups file (gf_nvt_pimd1), 
# and then 4 steps using the generated restrt files (nvt*.rst) as coordinate files in the groups file (gf_nvt_pimd2)

nc=24
# write gf_nvt_pimd1
for (( i=1; i<=${nc}; i++ )); do
    echo "-O -i pimd_nvt.in -p wat216.prmtop -c npt${i}.rst -r nvt${i}.rst" 
done > gf_nvt_pimd1

# run the first step
mpirun -np ${nc} sander.MPI -ng ${nc} -groupfile gf_nvt_pimd1
for (( i=1; i<=${nc}; i++ )); do
    cp nvt${i}.rst nvt${i}.rst.0 # backup the restrt files
done
cp pimdout pimdout.0 # backup pimdout file

# write gf_nvt_pimd2
for (( i=1; i<=${nc}; i++ )); do
    echo "-O -i pimd_nvt.in -p wat216.prmtop -c nvt${i}.rst -r nvt${i}.rst.new" 
done > gf_nvt_pimd2

nstep=4
for (( j=1; j<=${nstep}; j++ )); do
    # each run is a restart from previous step
    mpirun -np ${nc} sander.MPI -ng ${nc} -groupfile gf_nvt_pimd2
    for (( i=1; i<=${nc}; i++ )); do
        cp nvt${i}.rst.new nvt${i}.rst.${j} # backup the restrt files
        mv nvt${i}.rst.new nvt${i}.rst # prepare the restrt files for next run        
    done
    cp pimdout pimdout.${j}  # backup pimdout file
done
