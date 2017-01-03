;;; the distributed REPL

(in-package :mpi-worker)

(defun sequential-write (filename format &rest args)
  (loop for rank below (mpi-comm-size) do
    (when (= rank (mpi-comm-rank))
      (with-open-file (file filename :direction :output :if-exists :append)
        (apply #'format file format args)
        (force-output file)))
    (mpi-barrier)))

(defun main ()
  (mpi-init)
  (pushnew
   (lambda ()
     (mpi-finalize)
     (uiop:quit 0))
   sb-ext:*exit-hooks*)
  (loop for rank below (mpi-comm-size) do
    (when (= rank (mpi-comm-rank))
      (format t "~A~%" (mpi-comm-rank))
      (force-output)))
  (force-output)
  (let ((host (uiop/os:hostname))
        (port (+ 20000 (mpi-comm-rank)))
        (filename "connect-to-mpi-workers.lisp"))
    (when (= 0 (mpi-comm-rank))
      (with-open-file (file filename :direction :output :if-exists :supersede)
        (format file "(setf lfarm:*kernel*~%  (lfarm:make-kernel~%'(")
        (force-output file)))
    (mpi-barrier)
    (sequential-write filename "~%(~S ~A)" host port)
    (when (= 0 (mpi-comm-rank))
      (with-open-file (file filename :direction :output :if-exists :append)
        (format file ")))")
        (force-output file)))
    (mpi-barrier)
    (let ((pos (position "--script" sb-ext:*posix-argv* :test #'string-equal)))
      (if pos
          (load (elt sb-ext:*posix-argv* (1+ pos)))
          (lfarm-server:start-server host port)))))

