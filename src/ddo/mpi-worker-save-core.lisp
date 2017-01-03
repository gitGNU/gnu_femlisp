(in-package :mpi-worker)

(sb-ext:save-lisp-and-die
 "mpi-worker"
 :executable t
 :toplevel #'main)
