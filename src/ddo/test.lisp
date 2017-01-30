(in-package :ddo-test)

;;; This file assumes that you have compiled a binary mpi-worker
;;; in the directory #p"femlisp:bin" and started it with
;;; @code{mpirun -np xx mpi-worker}, where xx denotes the number
;;; of processes.  From a separate CL with DDO loaded you may then
;;; interactively evaluate the commands in this file.

(lfarm:end-kernel)
(worker-connect #p"femlisp:bin;mpi-worker-connection-data")

(ddo (fl.port:dynamic-space-size))

(ddo (in-package :ddo-test))
(ddo *package*)
(ddo (symbol-package (find-symbol "DUMMY")))

;;; number of workers
(defparameter *n* (lfarm:kernel-worker-count))
(ddo (defparameter *n* (mpi-comm-size)))
(ddo *n*)

(ddo *n*)
(ddo (synchronize))
(ddo (dbg-on :distribute))
(ddo (dbg-on :communication))
(ddo (dbg-on :local-synchronize))
(ddo (setq ddo::*debug-show-data* :all))
(ddo (distributed-data))
(ddo (synchronize))
(ddo (reset-distributed-objects))
(ddo (sb-ext:gc :full t))

;; an experiment with data merging
(ddo
  (defclass data ()
    ((data :initarg :data))))

(ddo (defparameter *x*
       (make-distributed-object
        (make-instance 'data :data (make-double-float-array 1 10.0))
        '(0 1)
        '(data))))

(ddo (distributed-slots *x*))
(ddo (insert-into-changed *x*))
(ddo (slot-value *x* 'data))
(ddo (synchronize))
(ddo (slot-value *x* 'data))

;; now, initialize *x* again, and change it on worker 1
(ddo (when (= (mpi-comm-rank) 1)
       (aref (slot-value *x* 'data) 0)))
(ddo (when (= (mpi-comm-rank) 1)
       (setf (aref (slot-value *x* 'data) 0) 3.0d0)
       (insert-into-changed *x*)))
(ddo (dbg-on :local-synchronize))
(ddo (let ((*synchronization-merger* #'minimum-id-merger))
       (synchronize)))
(ddo (slot-value *x* 'data))

;; some work on more  objects
(ddo (defparameter *x*
       (loop for k below 100
             collect
             (make-distributed-object
              (make-instance 'data :data (vector (coerce k 'double-float)))
              '(0 1))))
     (synchronize))
;; forget some objects on some workers
(ddo (setq *x* (remove-if (lambda (x)
                            (= (mod (aref (slot-value x 'data) 0)
                                    (mpi-comm-size))
                               (mpi-comm-rank)))
                          *x*)))
;; check lists
(ddo (length *x*))
;; and forget all
(ddo (setq *x* nil)
     (sb-ext:gc :full t)
     (synchronize))
(ddo (distributed-data))
