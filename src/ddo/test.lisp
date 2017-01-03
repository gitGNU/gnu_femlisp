(in-package :ddo-test)

;; (lfarm:end-kernel)

(load "connect-to-mpi-workers.lisp")

;; !!!  important: evaluate also the following:
(ddo (in-package :ddo-test))
(ddo *package*)
(ddo (symbol-package (find-symbol "DUMMY")))

;;; number of workers
(defparameter *n* (lfarm:kernel-worker-count))
(ddo (defparameter *n* (mpi-comm-size)))

;;; useful commands

(defun testing ()
  (ddo *n*)
  (ddo (unintern 'ddo-test::insert-into-changed))
  (ddo (load "packages.lisp"))
  (ddo- (load "ddo.lisp"))
  (ddo (synchronize))
  (ddo (accessing-exclusively ((count ddo::*local-id-count*))
         count))
  (ddo (accessing-exclusively ((it ddo::*new-distributed-objects*))
         (length it)))
  (ddo (make-instance 'ddo::dummy))
  (ddo (dbg-on :distribute))
  (ddo (dbg-on :communication))
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
          (make-instance 'data :data (fl.matlisp:double-vec 10.0))
          '(0 1)
          '(data))))
  
  (ddo (distributed-slots *x*))
  (ddo (insert-into-changed *x*))
  (ddo (slot-value *x* 'data))
  (ddo (synchronize))
  (ddo (slot-value *x* 'data))

  ;; the same using a container for the distributed objects which thus
  ;; share common mergers
  (ddo
    (defparameter *conti*
      (make-instance 'ddo-container-mixin
                     :distributed-slots (list (cons 'data (op-merger '+ 0))))))
  (ddo (defparameter *x*
         (make-distributed-object (make-instance 'data :data 10)
                                  '(0 1) *conti*)))
  (ddo (distributed-slots *x*))
  (ddo (insert-into-changed *x*))
  (ddo (slot-value *x* 'data))
  (ddo (synchronize))
  (ddo (slot-value *x* 'data))
  
  ;; then step through the following interaction:
  (ddo (distributed-data))
  (ddo (setq *x* nil))
  (ddo (distributed-data))
  (ddo (sb-ext:gc :full t))
  (ddo (distributed-data))
  (ddo (synchronize))
  (ddo (distributed-data))

  ;; now, initialize *x* again, and change it on worker 1
  (ddo (when (= (mpi-comm-rank) 1)
         (setf (slot-value *x* 'data) (vector 3))
         (insert-into-changed *x*)))
  ;; follow this with a (synchronize) and 
  (ddo (slot-value *x* 'data))
  
  ;; some work on more  objects
  (ddo (defparameter *x*
         (loop for k below 100
               collect
               (make-distributed-object (make-instance 'dummy :data (vector k))
                                        '(0 1 2 3))))
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
  )
