;;; -*- mode: lisp; fill-column: 70; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; synchronize.lisp - DDO synchronization
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

(file-documentation "DDO synchronization in two passes where the
first establishes the communication pattern and the second one
transfers the data.")

(defun all-processors ()
  (loop for k below (mpi-comm-size) collect k))

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
  "Default and simplest: all processors without me"
  (remove (mpi-comm-rank) (range< 0 (mpi-comm-size))))

(defmacro do-neighbors ((proc) &body body)
  `(loop for ,proc in (neighbors)
         do ,@body))

(defun first-comm-encode (vecs)
  "Encodes a list of integer vectors into a single UINT-VEC"
  (let ((nr-vectors (length vecs))
        (nr-values (reduce #'+ vecs :key #'length))
        (etype '(unsigned-byte 64)))
    (lret ((result (static-vectors:make-static-vector
                    (+ 1 nr-vectors nr-values)
                    :element-type etype
                    :initial-element 0)))
      (let ((k 0))
        (setf (aref result k) nr-vectors)
        (incf k)
        (loop for vec in vecs do
          (let ((vlength (length vec)))
            (setf (aref result k) (coerce vlength etype))
            (incf k)
            (inject! vec result k)
            (incf k vlength)))))))

(defun first-comm-decode (vec)
  "Decodes a vector encoded by FIRST-COMM-ENCODE into a list of integer vectors."
  (let ((len (length vec)))
    (and (plusp len)
         (let ((nr-vectors (aref vec 0))
               (k 1))
           (prog1
               (loop repeat nr-vectors
                 collect
                 (let ((vlength (aref vec k)))
                   (incf k)
                   (prog1
                       (subseq vec k (+ k vlength))
                     (incf k vlength))))
             (assert (= k len)))))))

(defun first-communication (new-objects changed-objects deleted-ids distribution)
  "A first communication sweep updates all distribution tables."
  (mpi-dbg :distribute "Entering first communication")
  (flet ((vec-of-lists ()
           (make-array (mpi-comm-size) :initial-element nil))
         (distribute (id place processors)
               (loop for proc in processors do
                 (unless (= proc (mpi-comm-rank))
                   (push id (aref place proc))))))
    (let ((new-to (vec-of-lists)) (changed-to (vec-of-lists)) (deleted-to (vec-of-lists))
          (first-comm-to (vec-of-lists)) (first-comm-from (vec-of-lists)))
      ;; new objects
      (loop for (object . processors) in new-objects do
        (distribute (local-id object) new-to processors))
      ;; changed objects
      (dotree (object changed-objects)
        (distribute (local-id object) changed-to (neighbors-for object)))
      ;; ids of deleted (garbage collected!) objects
      (dolist (id deleted-ids)
        (distribute id deleted-to (neighbors-for id)))
      ;; communication is to/from neighbors with whom old objects are shared
      ;; and for which new shared objects are generated
      ;; empty communication is therefore possible/necessary
      ;; and communication occuring is indicated by
      ;; (aref first-comm-for proc) being non-nil
      (do-neighbors (proc)
        (unless (= proc (mpi-comm-rank))
          (let ((new (aref new-to proc))
                (changed (aref changed-to proc))
                (deleted (aref deleted-to proc)))
          (when (or new changed deleted
                    (R-some distribution '_ proc '_)
                    *communicate-with-all*)
            (setf (aref first-comm-to proc)
                  (mapcar (rcurry #'coerce +ulong-vec+)
                          (list new changed deleted)))))))
      (assert (not (aref first-comm-to (mpi-comm-rank)))
              () "Should not communicate with myself")
      (mpi-dbg :distribute "Starting first communication")
      ;; parallel communication
      (let* ((comm-data
             (loop for proc in (neighbors)
                   when (aref first-comm-to proc)
                     collect (list :send proc (aref first-comm-to proc))
                     and collect (list :receive proc)))
           (comm-result (exchange comm-data :tag 1
                                            :encode #'first-comm-encode
                                            :cleanup #'static-vectors:free-static-vector
                                            :decode #'first-comm-decode)))
      (mpi-dbg (and *debug-show-data* :distribute) "Received ~A" comm-result)
      (loop for (from nil object) in comm-result
            do (setf (aref first-comm-from from)
                     object)))
    (mpi-dbg :distribute "Finished first communication")
    ;; all communication was successful, we can handle it:
    ;;
    ;; change distribution topology by inserting new and removing deleted objects
    (do-neighbors (proc)
      (awhen (aref first-comm-from proc)
        (destructuring-bind (distant-new distant-changed distant-deleted) it
          ;; distant-changed can be ignored (it is only relevant for the second step)
          (declare (ignore distant-changed))
          ;; insert new
          (let ((new (aref new-to proc)))
            (assert (= (length new) (length distant-new)))
            (loop for local-id in new
                  and distant-id across distant-new do
                    (R-insert distribution local-id proc distant-id)))
          ;; drop deleted
          (loop for distant-id across distant-deleted do
            (R-remove distribution '_ proc distant-id))
          )))
      ;; we pass all messages obtained in this first communication, because
      ;; we need part of it for interpreting the second communication
      first-comm-from)))

(deftype double-vec ()
  "Uniform @type{double-float} vector."
  '(simple-array double-float (*)))

(defun check-restricted-form (data)
  "Checks if data is of the form ((id double-vec ...) ...)."
  (let ((nr-objects 0)
        (nr-vectors 0)
        (nr-values 0))
    (loop for (id . vecs) in data
          do (incf nr-objects)
             (loop for vec in vecs do
               (incf nr-vectors)
               (assert (typep vec 'double-vec))
               (incf nr-values (length vec))))
    (values nr-objects nr-vectors nr-values)))

(defun inject! (x y offset)
  (loop for xc across x and i from offset do
    (setf (aref y i) xc)))

(defun second-comm-encode (data)
  "Encodes data of the form ((id double-vec ...) ...) into a double-vec."
  (multiple-value-bind (nr-objects nr-vectors nr-values)
      (check-restricted-form data)
    (lret ((result
             (static-vectors:make-static-vector
              (+ (* 2 nr-objects) nr-vectors nr-values)
              :element-type 'double-float
              :initial-element 0.0d0)))
      (let ((k 0))
        (loop for (id . vecs) in data do
          (setf (aref result k) (coerce id 'double-float))
          (incf k)
          (setf (aref result k) (coerce (length vecs) 'double-float))
          (incf k)
          (loop for vec in vecs do
            (let ((vlength (length vec)))
              (setf (aref result k) (coerce vlength 'double-float))
              (incf k)
              (inject! vec result k)
              (incf k vlength))))))))

(defun second-comm-decode (vec)
  "Decodes a vector encoded by SECOND-COMM-ENCODE."
  (loop with k = 0 and len = (length vec)
        until (>= k len)
        collect
        (multiple-value-bind (id frac) (floor (aref vec k))
          (assert (zerop frac))
          (incf k)
          (cons id
                (multiple-value-bind (nr-vectors frac) (floor (aref vec k))
                  (assert (zerop frac))
                  (incf k)
                  (loop repeat nr-vectors
                        collect
                        (multiple-value-bind (vlength frac) (floor (aref vec k))
                          (assert (zerop frac))
                          (incf k)
                          (prog1
                              (subseq vec k (+ k vlength))
                            (incf k vlength)))))))))

;;; (second-comm-decode (second-comm-encode (list (list 1 #d(2.0 3.0)) (list 2) (list 3 #d(2.0 3.0) #d()) )))

(defvar *fast-second-comm* t
  "If this variable is T, we introduce special encodings which avoid the
conspack library, but only work for slots with double vector data.")

(defun second-communication (new-objects changed-objects distribution first-comm-from)
  "The second communication step checks/modifies/unifies the slots of new and changed objects."
  (mpi-dbg :distribute "Entering second communication")
  (let ((second-comm-to (make-array (mpi-comm-size) :initial-element nil))
        (second-comm-from (make-array (mpi-comm-size) :initial-element nil)))
    ;; generate second messages
    (flet ((push-data (object)
             (let ((datum (cons (local-id object)
                                (distributed-slot-values object))))
               (loop for proc in (neighbors-for object) do
                 (push datum (aref second-comm-to proc))))))
      (dotree (object changed-objects)
        (push-data object))
      (dolist (object&processors new-objects)
        (push-data (car object&processors))))
    ;; communicate
    (let ((comm-result
            (apply #'exchange
                   (loop for proc in (neighbors)
                         when (aand (aref first-comm-from proc)
                                    (destructuring-bind (new changed deleted) it
                                      (declare (ignore deleted))
                                      (or (plusp (length new))
                                          (plusp (length changed)))))
                           collect (list :receive proc)
                         when (aref second-comm-to proc)
                           collect (list :send proc
                                         (aref second-comm-to proc)))
                   :tag 2
                   (when *fast-second-comm*
                     (list
                      :encode #'second-comm-encode
                      :cleanup #'static-vectors:free-static-vector
                      :decode #'second-comm-decode)))))
      (mpi-dbg (and *debug-show-data* :distribute) "second-comm-from=~A"
               comm-result)
      (loop for (from nil object) in comm-result
            do (setf (aref second-comm-from from)
                     ;; here: (decode-second-comm from object)
                     object)))
    ;; handle second receive
    (let ((change-table (make-hash-table)))
      ;; collect possible own changes into change-table
      (dotree (object changed-objects)
        (push (cons (mpi-comm-rank)
                    (distributed-slot-values object))
              (gethash (local-id object) change-table ())))
      (do-neighbors (proc)
        (let ((comm-from (aref second-comm-from proc)))
          (loop for (distant-id . slot-values) in comm-from
                for local-id = (caar (R-select distribution '_ proc distant-id))
                and k from 0 do
                  (assert local-id () "An entry (local-id ~D ~D) should exist on proc=~D!"
                          proc distant-id (mpi-comm-rank))
                  (let ((object (find-distributed-object local-id)))
                    (assert object () "Distributed object for local-id=~D should exist!" local-id)
                    ;; collect changed values
                    (push (cons proc slot-values)
                          (gethash local-id change-table ()))))))
      (when *synchronization-merger*
        (maphash (lambda (local-id changes)
                   (let ((object (find-distributed-object local-id)))
                     (mpi-dbg :local-synchronize "object=~A changes=~A" object changes)
                     (funcall *synchronization-merger* object changes)))
                 change-table)))
    ))

(defun synchronize (&optional label &rest args)
  "Synchronize distributed objects across the MPI kernel."
  ;; withouth MPI, synchronize does nothing
  (track-time (*synchronization-real-time*)
    (let ((times-string
            (format nil "~~&Synchronization~@[ in ~A~]~@[(~{~A~})~] needs ~~F seconds~%"
                    label args)))
      (measure-time-for-block (times-string :active-p (fl.debug:dbg-p :synchronize-timing))
        (when (mpi-initialized)
          ;; we establish exclusive access to all tables
          ;; at the moment we do not expect synchronization to be done concurrently with
          ;; other work, but later on this might be the case
          (accessing-exclusively ((new-objects *new-distributed-objects*)
                                  (changed-objects *changed-distributed-objects*)
                                  (deleted-ids *deleted-local-ids*)
                                  (distribution *distribution*))
            ;;
            ;; first communication
            ;;
            (let ((first-comm-from
                    (first-communication new-objects changed-objects deleted-ids
                                         distribution)))
              ;; the deleted local-ids have been communicated and should be handled
              ;; in the neighbors, so we clear the distribution data here, too
              (loop for id = (pop deleted-ids)
                    while id do (R-remove distribution id '_ '_))
              (mpi-dbg :distribute "Terminated first communication")
      
              (second-communication new-objects changed-objects distribution
                                    first-comm-from))
            (mpi-dbg :distribute "Terminated second communication")
    
            ;; finally we can clear also the new/changed objects tables
            (setf new-objects ()
                  changed-objects (make-binary-tree :red-black #'< :key #'local-id))
            (values)))))))
