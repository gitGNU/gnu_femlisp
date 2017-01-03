(in-package :mpi-extensions)

(defun random-double-vec (n)
  (let ((result
         (make-array n
                     :element-type 'double-float
                     :initial-element 0.0d0)))
    (loop for k below n do
         (setf (aref result k) (random 1.0d0)))
    result))

(defun random-messages (level nr-procs)
  (let ((n (expt 2 level)))
    (loop repeat nr-procs collect
         (random-double-vec n))))

(defun copy-to-static-vector (vec)
  (static-vectors:make-static-vector
   (length vec)
   :element-type (upgraded-array-element-type (array-element-type vec))
   :initial-contents vec))

(defmacro do-sequentially (&body body)
  (let ((rank (gensym "RANK")))
    `(loop for ,rank below (mpi-comm-size) do
      (when (= ,rank (mpi-comm-rank))
        ,@body)
      (mpi-barrier))))

(defun test (level repetitions)
  (do-sequentially
    (room)
    (sb-ext:gc :full t)
    (room)
    (format t "Proc: ~D, Level: ~D~%" (mpi-comm-rank) level)
    )
  (loop for k below repetitions
     do
       (let* ((messages (random-messages level (mpi-comm-size)))
              (comm-result
               (apply #'mpi-extensions:mpi-waitall-anything
                      (loop for proc below (mpi-comm-size)
                         for message = (elt messages proc)
                         unless (= proc (mpi-comm-rank))
                         collect
                           (mpi-extensions:mpi-isend-anything
                            message proc :tag 1
                            :encode #'copy-to-static-vector
                            :cleanup #'static-vectors:free-static-vector)
                         and collect
                           (mpi-extensions:mpi-irecv-anything
                            proc :tag 1 :decode #'copy-seq)))))
         (loop for (from nil object) in comm-result do
              (unless (equalp object (elt messages (mpi-comm-rank)))
                (format t "Proc ~D: Problem at level=~D&iteration=~D: received bad result from proc ~D~%"
                        (mpi-comm-rank) level k from)
                (error "Stop"))))))

(loop for level below 20 do
     (test level 1000))

(format t "&Proc ~D: Finishing program~%" (mpi-comm-rank))
