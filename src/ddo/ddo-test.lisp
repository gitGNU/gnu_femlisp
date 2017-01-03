(in-package :ddo-femlisp-test)

(defun random-double-vec (n)
  (lret ((result (make-double-vec n)))
    (loop for k below n do
      (setf (aref result k) (random 1.0)))))

(defun random-messages (level nr-procs &optional (order 5) (ncomps 3))
  (let ((nf (expt 4 level))
        (ne (* 4 (expt 4 level)))
        (nv (* 4 (expt 4 level))))
    (let* ((s (- order 1))
           (uf (* s s ncomps))
           (ue (* s ncomps))
           (uv (* ncomps))
           (total (+ (* 3 nf) (* nf uf)
                     (* 3 ne) (* ne ue)
                     (* 3 nv) (* nv uv))))
      (assert (< total 10e6) () "Messages probably too large")
      (let ((k 0))
        (loop repeat nr-procs
              collect
              (append
               (loop repeat nv
                     collect (list (incf k) (random-double-vec uv)))
               (loop repeat ne
                     collect (list (incf k) (random-double-vec ue)))
               (loop repeat nf
                     collect (list (incf k) (random-double-vec uf)))))))))

(defmacro do-sequentially (&body body)
  (with-gensyms (rank)
    `(loop for ,rank below (mpi-comm-size) do
      (when (= ,rank (mpi-comm-rank))
        ,@body)
      (mpi-barrier))))

(defun test (level repetitions)
  (do-sequentially
    (room)
    (sb-ext:gc :full t)
    (room)
    (dbg :ddo-test "Proc: ~D, Level: ~D~%" (mpi-comm-rank) level)
    )
  (loop for k below repetitions do
       (let ((messages (random-messages level (mpi-comm-size))))
         (let ((comm-result
                (ddo::exchange
                 (loop for proc below (mpi-comm-size)
                    for message = (elt messages proc)
                    when (not (= proc (mpi-comm-rank)))
                    nconc (list (list :receive proc)
                                (list :send proc message)))
                 :tag 2
                 :encode #'ddo::second-comm-encode
                 :cleanup #'static-vectors:free-static-vector
                 :decode #'ddo::second-comm-decode)))
           (loop for (from nil object) in comm-result do
                (unless (equalp object (elt messages (mpi-comm-rank)))
                  (format t "Proc ~D: Problem at level=~D&iteration=~D: received bad result from proc ~D~%"
                          (mpi-comm-rank) level k from)
                  (error "Stop")))))))

(defun random-messages-2 (level nr-procs)
  (let ((n (* 3 (expt 4 (1+ level)) 5)))
    (loop repeat nr-procs collect
         (random-double-vec n))))

(defun test-2 (level repetitions)
  (do-sequentially
    (room)
    (sb-ext:gc :full t)
    (room)
    (dbg :ddo-test "Proc: ~D, Level: ~D~%" (mpi-comm-rank) level)
    )
  (ignore-errors
    ;; simply decoding/encoding
    (loop for k below repetitions do
         (let ((messages (random-messages-2 level (mpi-comm-size))))
           (loop for proc below (mpi-comm-size)
              for message = (elt messages proc)
              for encoded = (ddo::copy-to-static-vector message)
              for decoded = (copy-seq encoded)
              do
                (assert (equalp message decoded))
                (static-vectors:free-static-vector encoded))))
    ;; then with communication
    (loop for k below repetitions do
         (let ((messages (random-messages-2 level (mpi-comm-size))))
           (let ((comm-result
                  (ddo::exchange
                   (loop for proc below (mpi-comm-size)
                      for message = (elt messages proc)
                      when (not (= proc (mpi-comm-rank)))
                      nconc (list (list :receive proc)
                                  (list :send proc message)))
                   :tag 2
                   :encode #'ddo::copy-to-static-vector
                   :cleanup #'static-vectors:free-static-vector
                   :decode #'copy-seq)))
             (loop for (from nil object) in comm-result do
                  (unless (equalp object (elt messages (mpi-comm-rank)))
                    (format t "Proc ~D: Problem at level=~D&iteration=~D: received bad result from proc ~D~%"
                            (mpi-comm-rank) level k from)
                    (error "Stop"))))))
    t))

(defun random-messages-3 (level nr-procs)
  (let ((n (expt 2 level)))
    (loop repeat nr-procs collect
         (random-double-vec n))))

(defun test-3 (level repetitions)
  (do-sequentially
    (room)
    (sb-ext:gc :full t)
    (room)
    (dbg :ddo-test "Proc: ~D, Level: ~D~%" (mpi-comm-rank) level)
    )
  (ignore-errors
    (loop for k below repetitions do
         (let ((messages (random-messages-3 level (mpi-comm-size))))
           (let ((comm-result
                  (ddo::exchange
                   (loop for proc below (mpi-comm-size)
                      for message = (elt messages proc)
                      when (not (= proc (mpi-comm-rank)))
                      nconc (list (list :receive proc)
                                  (list :send proc message)))
                   :tag 2
                   :encode #'ddo::copy-to-static-vector
                   :cleanup #'static-vectors:free-static-vector
                   :decode #'copy-seq)))
             (loop for (from nil object) in comm-result do
                  (unless (equalp object (elt messages (mpi-comm-rank)))
                    (format t "Proc ~D: Problem at level=~D&iteration=~D: received bad result from proc ~D~%"
                            (mpi-comm-rank) level k from)
                    (error "Stop"))))))
    t))

(defun test-4 (level repetitions)
  (do-sequentially
    (room)
    (sb-ext:gc :full t)
    (room)
    (dbg :ddo-test "Proc: ~D, Level: ~D~%" (mpi-comm-rank) level)
    )
  (ignore-errors
    (loop for k below repetitions
       do
         (let* ((messages (random-messages-3 level (mpi-comm-size)))
                (comm-result
                 (apply #'mpi-extensions::mpi-waitall-anything
                        (loop for proc below (mpi-comm-size)
                           for message = (elt messages proc)
                           unless (= proc (mpi-comm-rank))
                           collect
                             (mpi-extensions:mpi-isend-anything
                              message proc :tag 1
                              :encode #'ddo::copy-to-static-vector
                              :cleanup #'static-vectors:free-static-vector)
                           and collect
                             (mpi-extensions:mpi-irecv-anything
                              proc :tag 1 :decode #'copy-seq)))))
           (loop for (from nil object) in comm-result do
                (unless (equalp object (elt messages (mpi-comm-rank)))
                  (format t "Proc ~D: Problem at level=~D&iteration=~D: received bad result from proc ~D~%"
                          (mpi-comm-rank) level k from)
                  (error "Stop")))))
    t))


;;; (lfarm:end-kernel)
;;; (load "../connect-to-mpi-workers.lisp")
;;; (ddo- (load "ddo-femlisp/ddo-test.lisp"))
;;; (ddo (dbg-on :ddo-test))
;;; (ddo (setf ddo::*debug-show-data* nil))
;;; (ddo (loop for level below 6 do (test level 1000)))
