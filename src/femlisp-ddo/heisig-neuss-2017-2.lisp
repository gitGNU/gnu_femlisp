;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; heisig-neuss-2017-2.lisp - Calculations for [MHeisig-NNeuss-2017]-II
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2017
;;; Nicolas Neuss, Friedrich-Alexander-Universitaet Erlangen-Nuernberg
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions
;;; are met:
;;; 
;;; 1. Redistributions of source code must retain the above copyright
;;; notice, this list of conditions and the following disclaimer.
;;; 
;;; 2. Redistributions in binary form must reproduce the above copyright
;;; notice, this list of conditions and the following disclaimer in the
;;; documentation and/or other materials provided with the distribution.
;;; 
;;; THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED
;;; WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
;;; IN NO EVENT SHALL THE AUTHOR, THE FAU ERLANGEN-NUERNBERG, OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :cl-user)

(defpackage "FL.KONWIHR-PARALLEL"
  (:use "COMMON-LISP"
	"FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.DEBUG" "FL.TESTS" "FL.DEMO" "FL.PARALLEL"
        "FL.DICTIONARY"        
	"FL.FUNCTION" "FL.MESH"
	"FL.PROBLEM" "FL.CDR" "FL.ELLSYS" "FL.ELASTICITY" "FL.NAVIER-STOKES"
	"FL.DISCRETIZATION" "FL.ELLSYS-FE" "FL.ELASTICITY-FE" "FL.NAVIER-STOKES-FE"
	"FL.ITERATION" "FL.MULTIGRID" "FL.GEOMG"
	"FL.STRATEGY" "FL.PLOT"
	"FL.DOMAINS" "FL.APPLICATION" "FL.KONWIHR"
        "LFARM" "MPI"
        "RELATIONS" "NET.SCIPOLIS.GRAPHS" "DDO" "FEMLISP-DDO"
        "MPI-WORKER")
  (:export)
  (:documentation "Special package for the KONWIHR demo problem."))

(in-package :fl.konwihr-parallel)

(file-documentation
 "Demos for the MPI calculations reported in the KONWIHR paper by M. Heisig and N. Neuss.")

(in-package :fl.konwihr-parallel)

(defvar *mpi-workers* nil)
(defvar *mpi-konwihr-initialized* nil)
(defvar *mpi-konwihr-speed* nil)

(defun connect-to-mpi-workers
    (&optional (connection-spec #p"femlisp:bin;mpi-worker-connection-data"))
  (lfarm:end-kernel)
  (setq *mpi-workers* (worker-connect connection-spec))
  )

(defun disconnect-from-mpi-workers ()
  ;;(ddo (fl.port::quit))
  (lfarm:end-kernel)
  (setq *mpi-workers* nil)
  )

;;; (disconnect-from-mpi-workers)
;;; (connect-to-mpi-workers)

(defun konwihr-paper-max-levels-mpi ()
  (unless *mpi-workers*
    (error "No MPI workers"))
  (let ((mem (* (reduce #'min (ddo (fl.port:dynamic-space-size)))
                (lfarm:kernel-worker-count))))
    (floor (+ 1.9 (log (/ mem 1e9) 8)))))

;;; (konwihr-paper-max-levels)

(defun mpi-worker-connect-demo ()
  (when *mpi-workers*
    (format t "We disconnect from the existing pool of MPI workers...")
    (disconnect-from-mpi-workers)
    (setq *mpi-workers* nil
          *mpi-konwihr-speed* nil
          *mpi-konwihr-initialized* nil)
    (format t " Done.~%~%")
    )
  (format t "Please enter the connection data.  If the MPI worker pool
has been started on your local computer in the Femlisp
directory by 'make mpirun NP=xx', simply pressing ENTER might
already work.  Otherwise, you should enter either a string denoting a
filename on your local computer or the connection data by
copy-and-pasting the content of the file 'mpi-worker-connection-data'
written by your MPI workers at startup on the remote host.~%~%")
  (flet ((check-text (text)
           (or (and (zerop (length text)) :empty)
               (let* (*read-eval*
                      (object (read-from-string text)))
                 (cond ((listp object) :connection-data)
                       ((stringp object)
                        (unless (probe-file object)
                          (format t "File ~A not found.~%" text))
                        :filename)
                       (t (format t "Illegal input: ~A.~%" text)))))))
    (let ((connection-text (user-input-textfield #'check-text)))
      (connect-to-mpi-workers
       (case (check-text connection-text)
         (:empty nil)
         (:filename (probe-file connection-text))
         (:connection-data connection-text))))
    (cond
      (*mpi-workers*
       ;; we measure the speed of all workers
       (format t "~&Estimating CPU speed. This may take a little time, so please wait...")
       (force-output)
       (setq *mpi-konwihr-speed* (reduce #'min (ddo (konwihr-speed))))
       (format t "~&Measured ~D MFLOPS for slowest worker.~%"
               (* 10 (round *mpi-konwihr-speed* 10)))
       (force-output)
       (format t "We now initialize FEs and interpolation matrices for
the KONWIHR calculations on the workers.

This calculation may take about ~D seconds, but it is necessary only
once for all demos, so please wait..."
               (round (/ 20000 *mpi-konwihr-speed*)))
       (ddo (initialize-konwihr-paper-calculation))
       (format t "  Done.~%~%MPI workers are ready.~%"))
      (t
       (format t "~&Could not initialize MPI workers.~%")))
    (force-output t)))

(defun distributed-memory-demo ()
  (unless *mpi-workers*
    (format t "You must initialize an MPI worker pool before running this demo.")
    (return-from distributed-memory-demo))
  (let ((nr-workers (lfarm:kernel-worker-count))
        (allowed-worker-numbers '(2 3 4 6 8 12 16 24 48)))
    (unless (member nr-workers allowed-worker-numbers)
      (format t "For this demo we expect the number of MPI workers to
be an element of the set {~{~A~^,~}} for achieving an optimal load
distribution.  Please allocate a more appropriate worker kernel."
              allowed-worker-numbers)
      (return-from distributed-memory-demo)))
  (unless *mpi-konwihr-initialized*
    (format t "We initialize FEs and interpolation matrices
on the workers for the KONWIHR calculations, please wait...")
    (ddo (initialize-konwihr-paper-calculation))
    (setq *mpi-konwihr-initialized* t)
    (format t "  Done.~%~%"))
  ;; data initialization
  ;; GC
  (format t "For having more reproducible timing results, we do a
manual GC on all workers and clear all distributed data ...")
  (ddo
    (setq fl.application::*result* nil) ; drop previous results
    (sb-ext:gc :full t)
    (synchronize)
    (reset-distributed-objects))
  (format t "  Done.~%~%")
  (format t "Initialization phase complete - ready for calculating.~%~%")
  (let* ((nr-workers (lfarm:kernel-worker-count))
         (max-levels (konwihr-paper-max-levels-mpi))
         (query (format nil "Levels (1-~D): " (or max-levels 4)))
         (levels (user-input query #'parse-integer (_ (<= 1 _ max-levels))))
         (initial-mesh-refinements
           (if (member nr-workers '(2 3 6)) 0 1)))
    ;; stop multithreading on workers
    (ddo (end-kernel))
    ;; and allow it in some situations on request...
    (when (every (_ (> _ 1)) (ddo (length (fl.parallel::get-workers))))
      (when (eql (user-input
                  "Use threads on the workers (y/n)? "
                  (_ (cond
                       ((string-equal _ "y") :y)
                       ((string-equal _ "n") :n)))
                  (_ (member _ '(:y :n))))
                 :y)
        (ddo (progn (new-kernel) nil)  ; for dropping the return argument
             )))
    (format t "Sending commands to the workers and waiting for the result...")
    (force-output t)
    (let ((messages
            (eval
             `(ddo-capture
                (fl.parallel::end-kernel)
                ;; make sure that everything is initialized
                (initialize-konwihr-paper-calculation)
                ;; do the distributed calculation
                (let ((*distribute-n* ,nr-workers))
                  (heisig-neuss-2017-demo
                   (elasticity-inlay-cell-problem (n-cell-with-ball-hole 3))
                   :order 5 :levels ,levels :output 1 :distributed-p t
                   :initial-mesh-refinements ,initial-mesh-refinements))))))
      (format t "~&Got it!  Here is the worker's output:~%~A~%"
              (aref messages 0)))
    )
  )


(let ((worker-connect-demo
        (make-demo
         :name "Connect to workers"
         :short "Connect to a running worker pool"
         :long "When you have compiled and started a
Femlisp worker pool, you may connect to it using this demo.
Note that a local worker pool may be created and started by
issuing the commands

@code
make mpi-worker
make mpirun NP=<number-of-workers>
@end code

from inside the Femlisp directory."
         :execute (_ (mpi-worker-connect-demo))))
      (distributed-calculation-demo
        (make-demo
         :name "Distributed-memory parallel"
         :short "Run a distributed calculation"
         :long "Starts the calculation for calculating the
           effective elasticity tensor on the Femlisp worker
           pool."
         :execute (_ (distributed-memory-demo))
         :test-input "1~%"))
      (worker-disconnect-demo
        (make-demo
         :name "Disconnect from workers"
         :short "Disconnect from the worker pool"
         :long
         "When you are done with the tests, you may disconnect
from the worker pool using this 'demo'."
         :execute (_ (disconnect-from-mpi-workers)
                     (format t "Done")))))
  (adjoin-demo worker-connect-demo *heisig-neuss-2017-demo*)
  (adjoin-demo distributed-calculation-demo *heisig-neuss-2017-demo*)
  (adjoin-demo worker-disconnect-demo *heisig-neuss-2017-demo*)
  )

#|
(disconnect-from-mpi-workers)
(mpi-worker-connect-demo)
(distributed-memory-demo)
|#
