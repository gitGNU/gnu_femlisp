#!/bin/sh

echo "Zahl der Knoten:" $((2 * $PBS_NUM_NODES))

case `hostname | cut -b 1-1` in
    e) mpirun -np $((2 * $PBS_NUM_NODES)) -npernode 2 mpi-worker --dynamic-space-size 32000 --script ddo-femlisp/test-elahom-script.lisp
       ;;
    l) mpirun -np $((2 * $PBS_NUM_NODES)) -npernode 2 mpi-worker --dynamic-space-size 16000 --script ddo-femlisp/test-elahom-script.lisp
	;;
esac
