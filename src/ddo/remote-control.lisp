;;; -*- mode: lisp; fill-column: 70; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; remote-control.lisp - Remote control of a DDO worker pool
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2017-
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

(in-package :ddo)

(file-documentation
 "Useful functions for remote control of a worker pool.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  For interactive work
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *on-controller-p* nil
  "For testing purposes: If T perform ddo commands on the controller.")
  
(defmacro ddo (&body commands)
  "Broadcast the given commands"
  `(if *on-controller-p*
       (progn ,@commands)
       (lfarm:broadcast-task
        '(lambda () ,@commands))))

(defmacro ddox (&body commands)
  "Perform the commands serially on the controller"
  `(progn ,@commands))

(defmacro ddo- (&body commands)
  "Call workers sequentially, mainly for debugging purposes if separate
output is desired.  Of course, this should only be called with commands
which do not require synchronization!"
  `(if *on-controller-p*
       (progn ,@commands)
       (progn
         ,@(loop for rank below (lfarm:kernel-worker-count) collect
                 `(ddo
                    (when (= (mpi-comm-rank) ,rank)
                      (format t "~%>>> Working on processor ~D <<<" ,rank)
                      ,@commands
                      (force-output))))
         nil)))

