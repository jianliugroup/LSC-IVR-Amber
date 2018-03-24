#!/bin/bash
nsamp=3000
nc=24
for (( i=1; i<=${nsamp}; ++i ))
do
    mkdir ${i} # make work dir
    mpirun -n ${nc} sander.MPI -ng ${nc} -groupfile gf_sample # run pimd
    for (( j=1; j<=${nc}; ++j ))
    do
        cp nvt${j}.rst.new ${i}/        # copy generated coordinate files to work dir
        mv nvt${j}.rst.new nvt${j}.rst  # for next run
    done
done
