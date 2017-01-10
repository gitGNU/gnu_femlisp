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
                            :if-exists :append
                            :if-does-not-exist :create)
        (apply #'format file format args)
        (force-output file)))
    (mpi-barrier)))

(defun worker-connect (&optional (filespec-or-stream
                                  #p"femlisp:bin;mpi-worker-connection-data"))
  "Connect with the workers which are read from the given stream.
Each line has to contain a string (the hostname) and a number (the port)."
  (typecase filespec-or-stream
    (stream 
     (setq lfarm:*kernel*
           (lfarm:make-kernel
            (loop for line = (read-line filespec-or-stream nil) while line
                  collect
                  (let* ((items
                           (nth-value 1 (cl-ppcre:scan-to-strings
                                         "\"\([^\"]*)\"\\s*\([0-9]*\)" line))))
                    (unless items
                      (error "Bad connection data"))
                    (let ((host (aref items 0))
                          (port (parse-integer (aref items 1))))
                      (unless (and host (stringp host) port)
                        (error "Bad connection data"))
                      (list host port)))))))
    (t (with-open-file (stream filespec-or-stream)
         (worker-connect stream)))))

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
    ;; overwrite a possibly existing connection data
    (when (= 0 (mpi-comm-rank))
      (with-open-file (file filename :direction :output :if-exists :supersede)
        (force-output file)))
    ;; write new connection data
    (sequential-write filename "~S ~A~%" host port)
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

