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

(defun mpi-dbg (id format-string &rest args)
  (let ((rank (mpi-comm-rank)))
    (when (or (eq *report-ranks* T)
              (member rank *report-ranks*))
      (apply #'dbg id
             (format nil "Rank ~D:~A" rank format-string)
             args))))

(defmacro track-time ((var) &body body)
  (with-gensyms (timespan result)
    `(lret (,result)
       (let ((,timespan
               (measure-time (_ (setq ,result (progn ,@body))) 1 t)))
         (when ,var
           (incf ,var ,timespan))))))

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


