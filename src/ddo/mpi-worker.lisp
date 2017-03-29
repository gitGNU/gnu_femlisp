;;; -*- mode: lisp; fill-column: 70; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; mpi-worker.lisp - startup of DDO workers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2016
;;; Marco Heisig, Nicolas Neuss, FAU Erlangen-Nuernberg
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

(in-package :mpi-worker)

(defun sequential-write (filename format &rest args)
  (loop for rank below (mpi-comm-size) do
    (when (= rank (mpi-comm-rank))
      (with-open-file (file filename
                            :direction :output
                            :if-exists :append)
        (apply #'format file format args)
        (force-output file)))
    (mpi-barrier)))

(defun worker-connect (&optional connection-spec)
  "Connect with the workers which are derived from the given
connection-spec.  This can be empty, a string which is interpreted as
a filename, or a list where Each line has to contain a string (the
hostname) and a number (the port)."
  (ensure connection-spec #p"femlisp:bin;mpi-worker-connection-data")
  (typecase connection-spec
    (list
     (setq lfarm:*kernel* (lfarm:make-kernel connection-spec)))
    (stream
     (let (*read-eval*)
       (worker-connect (read connection-spec))))
    (string (worker-connect (make-string-input-stream connection-spec)))
    (t ;; should specify a pathname
     (with-open-file (stream connection-spec)
       (worker-connect stream)))))

;;; an interface suitable for the standard case when only
;;; a single MPI worker pool is used

(defvar *mpi-workers* nil
  "NIL if no standard pool of mpi-workers is active, otherwise the
corresponding lfarm kernel.")

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

(defun main ()
  (mpi-init)
  (pushnew
   (lambda ()
     (mpi-finalize)
     (uiop:quit 0))
   sb-ext:*exit-hooks*)
  ;; report workers ready on the console
  (loop for rank below (mpi-comm-size) do
    (when (= rank (mpi-comm-rank))
      (format t "Worker with MPI rank = ~D ready~%" rank)
      (force-output))
    (mpi-barrier))
  ;;
  (let ((host (uiop/os:hostname))
        (port (+ 20000 (mpi-comm-rank)))
        (filename "mpi-worker-connection-data"))
    ;; overwrite a possibly existing connection data file
    (when (= 0 (mpi-comm-rank))
      (with-open-file (file filename :direction :output :if-exists :supersede)
        (format file "(~%")
        (force-output file)))
    ;; write new connection data
    (sequential-write filename "(~S ~A)~%" host port)
    ;; terminate
    (when (= 0 (mpi-comm-rank))
      (with-open-file (file filename :direction :output :if-exists :append)
        (format file ")~%")
        (force-output file)))
    ;; check command line
    (if (find-if (lambda (arg)
                 (member arg '("--script" "--eval") :test #'string=))
                 sb-ext:*posix-argv*)
        ;; process scripts/commands from the command line
        (loop for args on sb-ext:*posix-argv* by #'rest
            for arg = (first args) do
              (cond ((string= arg "--script")
                     (awhen (second args) (load it)))
                    ((string= arg "--eval")
                     (awhen (second args)
                       (eval (read-from-string it))))))
        ;; start a server waiting for remote control
        (lfarm-server:start-server host port))))

