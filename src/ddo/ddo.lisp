;;; -*- mode: lisp; fill-column: 75; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ddo.lisp - dynamic distributed objects
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2015-
;;; Nicolas Neuss, FAU Erlangen-Nuernberg
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
;;; IN NO EVENT SHALL THE AUTHOR, THE KARLSRUHE INSTITUTE OF TECHNOLOGY,
;;; OR OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :ddo)

;;; local ids

(defvar *local-id-count* (mutex-wrap 0))

(defun new-local-id ()
  (accessing-exclusively ((count *local-id-count*))
    (incf count)))

(defvar *distributed-objects*
  (mutex-wrap (make-hash-table :weakness :value))
  "A weak hash table mapping local-id to object.")

(defun find-distributed-object (id)
  (accessing-exclusively ((table *distributed-objects*))
    (gethash id table)))

(defvar *distribution*
  (mutex-wrap (make-number-relation 3))
  "A relation local-id<->processor<->distant-id.")

(defun distribution-list ()
  "For testing purposes look at the distribution entries"
  (accessing-exclusively ((R *distribution*))
    (R-select R '_ '_ '_)))

(defvar *deleted-local-ids* (mutex-wrap ())
  "A list of local-ids whose objects are not needed anymore on this processor")

(defvar *new-distributed-objects*
  (mutex-wrap ())
  "Association list mapping new distributed objects to its owners.")

(defgeneric local-id (object)
  (:method (object)
      "Default method for non-distributed objects returns NIL."
    nil))

(defvar *changed-distributed-objects*
  (mutex-wrap (make-binary-tree :red-black #'< :key #'local-id))
  "Table containing changed distributed objects.")

(defun reset-distributed-objects ()
  "For testing purposes!  Resets all distributed-object tables."
  (accessing-exclusively ((objects *distributed-objects*))
    (setf objects (make-hash-table :weakness :value)))
  (accessing-exclusively ((d *distribution*))
    (setf d (make-number-relation 3)))
  (accessing-exclusively ((ids *deleted-local-ids*))
    (setf ids ()))
  (accessing-exclusively ((objects *new-distributed-objects*))
    (setf objects ()))
  (accessing-exclusively ((objects *changed-distributed-objects*))
    (setf objects (make-binary-tree :red-black #'< :key #'local-id)))
  (values))

(defun distributed-data ()
  "Only for debugging purposes: return all distributed data in a property list"
  (accessing-exclusively ((new-objects *new-distributed-objects*)
                          (changed-objects *changed-distributed-objects*)
                          (deleted-ids *deleted-local-ids*)
                          (distribution *distribution*)
                          (distributed-objects *distributed-objects*))
    (list :new-objects new-objects
          :changed-objects (net.scipolis.relations::tree-leaves changed-objects)
          :deleted-ids deleted-ids
          :distribution (R-select distribution '_ '_ '_)
          :distributed-objects (ALEXANDRIA.0.DEV:HASH-TABLE-ALIST distributed-objects))))

;;; distributed object and distributed container datatypes

(defclass ddo-container-mixin ()
  ((distributed-slots
    :initform ()
    :initarg :distributed-slots
    :documentation "This determines a list of slots for which changes are
propagated.  This is a list consisting of items that are either a symbol or
a pair of the form (symbol . combiner).  Here, the symbol is the name of a
slot and the combiner is a certain function which is called when combining
distributed values."))
  (:documentation "Mixin for a ddo container which may contain distributed
and non-distributed objects.  Not used at the moment!"))

#+(or)
(defun distributed-container-p (object)
  (typep object 'ddo-container-mixin))

(defclass ddo-mixin (property-mixin)
  ((local-id
    :reader local-id
    :initform (new-local-id)
    :documentation "The local id of this object")
   (distributed-slots :reader distributed-slots :initform () :initarg :container :type list)))

(defun distributed-p (object)
  (typep object 'ddo-mixin))

(defun distributed-finalizer (local-id)
  "Returns a function which is called when the distributed object with
local-id equal to @arg{local-id} is garbage-collected."
  (lambda ()
    (accessing-exclusively ((ids *deleted-local-ids*))
      (push local-id ids))))

#+(or)
(defgeneric container (obj)
  (:documentation "Returns the container of the ddo @arg{obj} or nil if there is none.")
  (:method ((obj ddo-mixin))
      (let ((c/s (slot-value obj 'container-or-slots)))
        (and (typep c/s 'ddo-container-mixin)
             c/s))))

#+(or)
(defgeneric distributed-slots (object)
  (:documentation "Returns a list of slot-names or slot-name/combiner pairs.")
  (:method (object)
      "Default: No data has to be synchronized."
    nil)
  (:method ((object ddo-mixin))
      (let ((c/s (slot-value object 'container-or-slots)))
        (if (listp c/s)
            c/s
            (slot-value c/s 'distributed-slots)))))

(defun distributed-slot-names (object)
  (mapcar (lambda (slot)
            (etypecase slot
              (symbol slot)
              (cons (car slot))))
          (distributed-slots object)))

(defgeneric distributed-slot-values (object)
  (:method (object)
      "By default return nothing for non-distributed objects."
    ())
  (:method ((object ddo-mixin))
      "By default return all distributed slot values for distributed objects."
    (mapcar (curry #'slot-value object)
            (distributed-slot-names object))))

(defun ensure-distributed-class (class &optional (type :object))
  "Generates for @arg{class} a distributed object/container variant."
  (let ((mixin (ecase type
                 (:object 'ddo-mixin)
                 #+(or)
                 (:container 'ddo-container-mixin))))
    (if (subtypep class mixin)
        class
        (let ((name (class-name class)))
          (fl.amop:find-programmatic-class
           (list class mixin)
           :class-name (intern (format nil "=~A=" name) (symbol-package name)))))))

(defun make-distributed-object (object processors &optional container-or-slots)
  "Make OBJECT into a distributed object belonging to PROCESSORS.
The change will become active only after the next synchronization!"
  (assert (not (typep object 'ddo-mixin)))
  (assert (member (mpi-comm-rank) processors))
  (assert (every (rcurry #'< (mpi-comm-size)) processors))
  (when (> (length processors) 1)
    ;; (assert (apply #'< processors)) ; not necessary
    (change-class object
                  (ensure-distributed-class (class-of object))
                  :container container-or-slots)
    (sb-ext:finalize object (distributed-finalizer (local-id object)))
    (accessing-exclusively ((objects *distributed-objects*)
                            (new-objects *new-distributed-objects*))
      ;; establish new objects as distributed objects
      (setf (gethash (local-id object) objects) object)
      (push (cons object processors) new-objects)))
  object)

(defun make-distributed-container (object &optional distributed-slots)
  "Turn the container @arg{object} into a distributed container."
  (assert (not (typep object 'ddo-container-mixin)))
  (change-class object
                (ensure-distributed-class (class-of object) :container)
                :distributed-slots distributed-slots))

;;; A simple slot merger
(defgeneric minimum-id-merger (object id-value-pairs)
  (:documentation "In principle, this function should look up up the
 minimal id and set all distributed slots accordingly.")
  (:method (object id-value-pairs)
      "Do nothing on non-distributed objects.")
  (:method ((object ddo-mixin) id-value-pairs)
      "Look up the minimal id in @arg{id-value-pairs} and set the
distributed slots accordingly."
    (when id-value-pairs
      (let* ((min-id (loop for entry in id-value-pairs minimize (car entry)))
             (values (cdr (assoc min-id id-value-pairs)))
             (dslots (distributed-slots object)))
        (assert (= (length dslots) (length values)))
        (loop for slot in (distributed-slots object)
              and value in values
              for slot-name = (if (symbolp slot) slot (car slot))
              do
                 (setf (slot-value object slot-name) value))))))

(defun op-merger (reduce-op initial-value)
  (lambda (object id-value-pairs)
    (declare (ignore object))
    (reduce reduce-op id-value-pairs
            :key #'cdr :initial-value initial-value)))

(defvar *synchronization-merger* nil
  "Merger for the synchronization which can be a generic function
  dispatching on the class of object.  An example is provided by the
  generic function #'minimum-id-merger.  Usually this dynamic variable will
  be bound around a call to @func{synchronize}.")

(defclass dummy ()
  ((data :initform nil :initarg :data))
  (:documentation "For testing only."))

(defmethod print-object ((d dummy) stream)
  (format stream "DUMMY{~A}" (slot-value d 'data)))

(defun insert-into-changed (object)
  (assert (typep object 'ddo-mixin))
  (accessing-exclusively ((changed *changed-distributed-objects*))
    (insert object changed)))

;;; One could think about turning this into an after method for
;;; (setf slot-value) like shown here:

#+(or)
(defmethod (setf slot-value-using-class) :after
    (new-value (class sb-mop:standard-class) (object ddo-mixin)
               (sloctd sb-mop:standard-effective-slot-definition))
  (insert-into-changed object))

(cl-store:defstore-cl-store (obj function stream)
  (cl-store::store-simple-string "<some-function>" stream))

(cl-store:defstore-cl-store (obj standard-object stream)
  (cl-store::store-simple-string (symbol-name (class-name (class-of obj))) stream))

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

;;; (ddo-performance-check)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-distributed ()
  (let ((test (make-number-relation 3)))
    (finalize test
              (lambda ()
                (format *trace-output* "~&Finalizing~%")
                (force-output *trace-output*)))
    (gc :full t))
  (let ((test (make-number-relation 3)))
    (finalize test
              (lambda ()
                (format *trace-output* "~&Finalizing~%")
                (force-output *trace-output*)))
    (accessing-exclusively ((table *distributed-objects*))
      (setf (gethash 1 table) test))
    (gc :full t)
    (accessing-exclusively ((table *distributed-objects*))
      (print (hash-table-count table)))
    ;; entry is apparently cleared only later when (gc :full t) is called separately
    )

  (let ((object (make-instance 'dummy :data 10)))
    (ensure-distributed-class (find-class 'dummy))
    (make-distributed-object object '(0))
    (finalize object (distributed-finalizer (local-id object)))
    (gc :full t)
    )
  (dbg-on :distribute)
  (reset-distributed-objects)

  (accessing-exclusively ((distribution *distribution*))
    (pprint-tree (cdr (first (net.scipolis.relations::indices distribution)))))
  
  )


