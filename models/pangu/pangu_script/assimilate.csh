#!/bin/csh
#BSUB -J filter
#BSUB -q mpi
#BSUB -n 24
#BSUB -W 5
#BSUB -o log.filter
#BSUB -e err.filter

mpirun ./filter
if ( -e ${RUN_DIR}/obs_seq.final )  touch ${RUN_DIR}/filter_done

exit 0

