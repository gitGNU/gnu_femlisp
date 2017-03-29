;;; -*- mode: lisp; fill-column: 70; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; utils.lisp - DDO utility functions
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

(defun mpi-rank ()
  (when (mpi-initialized)
    (mpi-comm-rank)))

(defun mpi-dbg (id format-string &rest args)
  (let ((rank (mpi-comm-rank)))
    (when (or (eq *report-ranks* T)
              (member rank *report-ranks*))
      (apply #'dbg id
             (format nil "Rank ~3D:~A" rank format-string)
             args))))

(defun mpi-debug-file (&optional rank)
  (parse-namestring
   (format nil "femlisp:data;mpilog-~3,'0D"
           (or rank (mpi-comm-rank)))))

(defmacro with-split-mpi-debug-output (&body body)
  (with-gensyms (stream)
    `(with-open-file (,stream (mpi-debug-file)
                              :direction :output
                              :if-exists :append
                              :if-does-not-exist :create)
       (let ((fl.debug::*debug-io* ,stream))
         (format fl.debug::*debug-io* "~&*** Rank ~3D ***~%" (mpi-comm-rank))
         (force-output fl.debug::*debug-io*)
         ,@body))))
  
(defmacro track-time ((var) &body body)
  (with-gensyms (timespan result)
    `(lret (,result)
       (let ((,timespan
               (measure-time (_ (setq ,result (progn ,@body))) 1 t)))
         (when ,var
           (incf ,var ,timespan))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  Helper functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun neighbors-for (object-or-local-id)
  "All processors for a distributed object"
  (let ((local-id (etypecase object-or-local-id
                    (ddo-mixin (local-id object-or-local-id))
                    (integer object-or-local-id))))
    (accessing-exclusively ((it *distribution*))
      (mapcar #'second (R-select it local-id '_ '_)))))

(defun owners (object-or-local-id)
  (cons (mpi-comm-rank) (neighbors-for object-or-local-id)))

(defun masterp (object-or-local-id)
  "The master of a distributed object is simply as the processor with
minimal rank.  Note that using this function somehow works against the
whole DDO concept, so that its use should be avoided whenever possible."
  (apply #'< (mpi-comm-rank) (neighbors-for object-or-local-id)))

(defun neighbors ()
  "Look with which processors we are connected."
  (let ((rank (mpi-comm-rank)))
    (accessing-exclusively ((distribution *distribution*))
      (loop for proc below (mpi-comm-size)
            unless (or (= rank (mpi-comm-rank))
                       (not (R-select distribution '_ proc '_)))
              collect proc))))

(defmacro do-neighbors ((proc) &body body)
  `(loop for ,proc in (neighbors)
         do ,@body))

(defmacro do-processors ((proc) &body body)
  `(loop for ,proc below (mpi-comm-size)
         do ,@body))

(defun all-processors ()
  (loop for k below (mpi-comm-size) collect k))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  For profiling relation handling
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ddo-performance-check ()
  "Checks how much time the relation handling alone, that is without any
communication, needs on the currently active distributed objects."
  (accessing-exclusively ((distribution *distribution*)
                          (distributed-objects *distributed-objects*))
    (format t "~&DDO performance check~%")
    (format t "~&Number of distributed-objects = ~D~%" (hash-table-count distributed-objects))
    (measure-time-for-block ("~&Selecting all processor interfaces: ~F~%")
      (do-neighbors (proc)
        (R-select distribution '_ proc '_)))
    (let (entries)
      (measure-time-for-block ("~&Selecting all distribution entries: ~F~%" :active-p t)
        (setq entries (R-select distribution '_ '_ '_)))
      (format t "~&#Number of distribution-entries = ~D~%" (length entries))
      (measure-time-for-block ("~&Selecting all distribution entries sequentially: ~F" :active-p t)
        (loop for entry in entries do
          (apply #'R-select distribution entry)))
      (measure-time-for-block ("~&Selecting local-id sequentially: ~F" :active-p t)
        (loop for (local-id proc distant-id) in entries do
          (R-select distribution '_ proc distant-id)))
      (measure-time-for-block ("~&Selecting distant-id sequentially: ~F" :active-p t)
        (loop for (local-id proc distant-id) in entries do
          (R-select distribution local-id proc '_))))
    (format t "~&END performance check~%")
    (force-output)))

#+(or)
(defun send (data proc)
  (error "Not used anymore")
  (assert (/= proc (mpi-comm-rank)) ()
          "Sending to myself: proc=~D" proc)
  (mpi-dbg :distribute "Sending ~A to ~A"
           (if *debug-show-data* data :some-data)
           proc)
  (when (mpi-initialized)
    (track-time (*communication-real-time*)
      (mpi-send-anything data proc))))

#+(or)
(defun receive (proc)
  (error "Not used anymore.")
  (assert (/= proc (mpi-comm-rank)) ()
          "Receiving from myself: proc=~D" proc)
  (track-time (*communication-real-time*)
    (cond ((mpi-initialized)
           (lret ((result (third (first (mpi-recv-anything proc)))))
             (mpi-dbg :distribute "Received ~A from ~A"
                      (if *debug-show-data* result :some-data)
                      proc)))
          (t (if (zerop proc)
                 (error "Called receive from myself.")
                 (error "No other processor"))))))

(defun total-size (item)
  (etypecase item
    (number 1)
    (array (array-total-size item))
    (list (loop for elem in item sum (total-size elem)))))

(defun copy-to-static-vector (vec)
  (static-vectors:make-static-vector
   (length vec)
   :element-type (upgraded-array-element-type (array-element-type vec))
   :initial-contents vec))

(defun exchange (request-list
                 &key (tag 0) encode cleanup decode)
  "Wait for all requests in request-list and return a list of receive
  requests.  A request is either of the form
  (:send proc data) or (:receive proc)."
  (let ((data-size
          (and *communication-size*
               (loop for request in request-list
                     sum (total-size (third request))))))
    (when data-size
      (incf *communication-size* data-size))
    (dbg-when :communication
      ;; generate exchange list
      (let* ((n (mpi-comm-size))
             (from (make-array n :initial-element 0))
             (to (make-array n :initial-element 0)))
        (loop for (type proc nil) in request-list do
          (ecase type
            (:send (incf (aref to proc)))
            (:receive (incf (aref from proc)))))
        (mpi-dbg :communication "Exchange pattern: SEND=~{~1D~} RECEIVE=~{~1D~}"
                 (coerce to 'list) (coerce from 'list))))
    (mpi-dbg (and *debug-show-data* :communication) "Exchanging ~A"
             (cond ((eq *debug-show-data* :all) request-list)
                   (t (if data-size
                          (format nil "~D numbers" data-size)
                          :some-data)))))
  (track-time (*communication-real-time*)
    (apply #'mpi-waitall-anything
           (loop for request in request-list
                 collect
                 (destructuring-bind (type proc &optional data) request
                   (ecase type
                     (:send
                      (unless data
                        (error "Processor ~D: No data field for ~A-operation to ~D"
                               (mpi-comm-rank) type proc))
                      (apply #'mpi-isend-anything data proc :tag tag
                             (and encode (list :encode encode :cleanup cleanup))))
                     (:receive
                      (apply #'mpi-irecv-anything proc :tag tag
                             (and decode (list :decode decode))))))))))

(defun unsymmetric-at (array i j)
  (/= (aref array i j)
      (aref array j i)))

(defun array-symmetric-p (array)
  (let ((dims (array-dimensions array)))
    (assert (samep dims) () "Array is non-quadratic")
    (let ((n (first dims)))
      (dotimes (i n)
        (dotimes (j i)
          (when (unsymmetric-at array i j)
            (return-from array-symmetric-p (values nil i j))))))
    t))

(defun connection-string (connections)
  "Prints a connection array in compact quadratic form."
  (with-output-to-string (stream)
    (loop for i below (array-dimension connections 0) do
      (format stream "~{~2A~}~%"
              (loop for j below (array-dimension connections 1)
                    collect
                    (ecase (aref connections i j)
                      (0 0)
                      (1 (if (= (aref connections j i) 1)
                             1
                             '*))))))))

(defun distribution-consistency-check ()
  "Does a communication between all processors and checks if the distribution fits.
The result should be a vector containing only zeros."
  (accessing-exclusively ((distribution *distribution*)
                          (distributed-objects *distributed-objects*))
    (let ((sharing (map +ulong-vec+
                        (lambda (proc)
                          (if (R-select distribution '_ proc '_) 1 0))
                        (all-processors))))
      (let* ((comm-data
               (loop for proc in (all-processors)
                     collect (list :send proc sharing)
                     collect (list :receive proc)))
             (comm-result (exchange comm-data :tag 0))
             (my-rank (mpi-comm-rank)))
        (let* ((items (cons (list my-rank nil sharing)
                            comm-result))
               (connections (make-array (twice (mpi-comm-size)))))
          ;; fill with content
          (loop for (proc nil data) in items do
            (do-processors (k)
              (setf (aref connections proc k)
                    (aref data k))))
          (let ((consistent-p (array-symmetric-p connections)))
            (unless consistent-p
              (mpi-dbg :communication "Inconsistent connections:~%~A"
                       (connection-string connections))
              ;; Report own inconsistencies if there should be any
              (do-processors (proc)
                (when (and (unsymmetric-at connections my-rank proc)
                           (eql (aref connections my-rank proc) 1))
                  (mpi-dbg :communication "~&For processor ~D:~%" proc)
                  (loop for (local-id nil nil)
                          in (R-select distribution '_ proc '_)
                        do
                           (mpi-dbg :communication "~A~%"
                                    (gethash local-id distributed-objects))))))
            (values consistent-p connections)))))))

;;; Debugging

(defun analyze-mpi-debug-output ()
  "Analyzes send/receive output.  This is probably done better prior
to synchronization by distribution-consistency-check."
  (let* ((np (lfarm:kernel-worker-count))
         (outputs
           (loop for rank below np collect
             (with-open-file (stream (ddo::mpi-debug-file rank))
               (loop for line = (read-line stream nil nil)
                     while line collect line)))))
    (assert (samep (mapcar #'length outputs)))
    (let ((nr-comms (length (first outputs))))
      (let ((send-pattern (make-array (list nr-comms np np) :initial-element nil))
            (receive-pattern (make-array (list nr-comms np np) :initial-element nil)))
        (loop for comm below nr-comms do
          (loop for output in outputs
                and rank from 0
                for line = (nth comm output)
                when (search "SEND=" line) do
                  (let ((sends (+ 5 (search "SEND=" line)))
                        (receives (+ 8 (search "RECEIVE=" line))))
                    (loop for k below np do
                      (when (eql (aref line (+ k sends)) #\1)
                        (setf (aref send-pattern comm rank k) T))
                      (when (eql (aref line (+ k receives)) #\1)
                        (setf (aref receive-pattern comm rank k) T))))))
        ;; check that they are transposed in the latter variables
        (loop with flag = nil
              for comm below nr-comms
              do
                (loop for i below np do
                  (loop for j below np do
                    (unless (eql (aref send-pattern comm i j)
                                 (aref receive-pattern comm j i))
                      (format t "Problem at comm=~D: (~D ~A-~A ~D)~%"
                              comm i
                              (if (aref send-pattern comm i j) #\> #\|)
                              (if (aref receive-pattern comm j i) #\> #\|)
                              j)
                      (setf flag t))))
                (when flag
                  (return (list (lambda (i j) (aref send-pattern comm i j))
                                (lambda (i j) (aref receive-pattern comm i j))))))))))


(defun test-utils ()
  (ddo-performance-check)
  (make-array '(2 2) :initial-contents '(#(1 2) #(3 4)))
  (connection-string #2a((0 1) (0 0)))
  (array-symmetric-p #2a((0 0) (0 0)))
  )

