(in-package :fl.parallel)

(error "This file should not be loaded.  It was once a quickhack
for working around a bad CPU affinity setting of mpirun.")

(defun set-cpu-affinity-to-mpi-comm-rank ()
  (when (cl-mpi:mpi-initialized)
    (let* ((rank (cl-mpi:mpi-comm-rank))
           (proc (nth rank (get-processors)))
           (id (pi-processor proc)))
      (format t "~&Setting affinity of rank ~D to ~D~%" rank id)
      (sb-cpu-affinity:with-cpu-affinity-mask (mask :save t)
        (sb-cpu-affinity:clear-cpu-affinity-mask mask)
        (setf (sb-cpu-affinity:cpu-affinity-p id mask) t)))))

(set-cpu-affinity-to-mpi-comm-rank)
