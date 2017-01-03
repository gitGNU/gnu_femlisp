#!/bin/sh

export TOTAL_MEMORY=240000
export NUM_SOCKS=4
export NUM_NODES=48
export N_PER_SOCK=$(($NUM_NODES / $NUM_SOCKS))
export N_THREADS=$((12 / $N_PER_SOCK))
export LOCAL_MEMORY=$(($TOTAL_MEMORY / $NUM_NODES))

echo "TOTAL_MEMORY=" $TOTAL_MEMORY
echo "NUM_SOCKS=" $NUM_SOCKS
echo "NUM_NODES=" $NUM_NODES
echo "N_PER_SOCK=" $N_PER_SOCK
echo "N_THREADS=" $N_THREADS
echo "LOCAL_MEMORY=" $LOCAL_MEMORY

#mpirun -np $NUM_NODES --map-by ppr:$N_PER_SOCK:socket:pe=$N_THREADS ./mpi-worker --dynamic-space-size $LOCAL_MEMORY --script ddo-femlisp/konwihr-paper.lisp
#mpirun -np $NUM_NODES --map-by ppr:$N_PER_SOCK:socket:pe=$N_THREADS ./mpi-worker --dynamic-space-size $LOCAL_MEMORY --script ddo-femlisp/test-elahom-script.lisp

mpirun -np 48 --map-by ppr:12:socket ./mpi-worker --dynamic-space-size 4000 --script ddo-femlisp/test-elahom-script.lisp
