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

(defparameter *debug-show-data* nil
  "When T show all data communicated between processors.  This is only reasonable for
toy problems.")

(defvar *communicate-with-all* t
  "For debugging purposes, this may be set to T, for ensuring that communication
is not blocked by a non-fitting communication pattern.")

(defvar *communication-real-time* nil
  "When non-nil, communication time is recorded.")

(defvar *synchronization-real-time* nil
  "When non-nil, communication time is recorded.")

(defvar *communication-size* nil
  "When non-nil, communication size is recorded.")

(defvar *report-ranks* t
  "List of ranks for which reports are shown.  T means all ranks.")

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

;;; Some slot mergers
(defun minimum-id-merger (object id-value-pairs)
  (declare (ignore object))
  (and id-value-pairs
       (let ((min-id (loop for entry in id-value-pairs minimize (car entry))))
         (cdr (assoc min-id id-value-pairs)))))

(defun op-merger (reduce-op initial-value)
  (lambda (object id-value-pairs)
    (declare (ignore object))
    (reduce reduce-op id-value-pairs
            :key #'cdr :initial-value initial-value)))

(defvar *synchronization-merger* nil
  "Merger for the synchronization which can be a generic function dispatching on the class of object, the object itself and the slots.")

(defclass dummy ()
  ((data :initform nil :initarg :data))
  (:documentation "For testing only."))

(defmethod print-object ((d dummy) stream)
  (format stream "DUMMY{~A}" (slot-value d 'data)))

(defun insert-into-changed (object)
  (assert (typep object 'ddo-mixin))
  (accessing-exclusively ((changed *changed-distributed-objects*))
    (insert object changed)))

;;; Turn into an after method for (setf slot-value)??

#+(or)
(defmethod (setf slot-value-using-class) :after
    (new-value (class sb-mop:standard-class) (object ddo-mixin)
               (slotd sb-mop:standard-effective-slot-definition))
  (insert-into-changed object))

;;; Flags coded as integers for a communication stream

(defconstant +deleted+ 0)
(defconstant +changed+ 1)
(defconstant +new+ 2)

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

(defvar +ulong+ '(unsigned-byte 64))
(defvar +ulong-vec+ `(simple-array ,+ulong+ (*)))

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
            (exchange
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
             :encode #'second-comm-encode
             :cleanup #'static-vectors:free-static-vector
             :decode #'second-comm-decode)))
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
                  ;;; (assert (= local-id distant-id))  ;; only for testing on 2 procs
                  (assert local-id () "An entry (local-id ~D ~D) should exist on proc=~D!"
                          proc distant-id (mpi-comm-rank))
                  (let ((object (find-distributed-object local-id)))
                    (assert object () "Distributed object for local-id=~D should exist!" local-id)
                    ;; collect changed values
                    (push (cons proc slot-values)
                          (gethash local-id change-table ()))))))
      (when *synchronization-merger*
        (maphash (lambda (local-id changes)
                   (funcall *synchronization-merger*
                            (find-distributed-object local-id)
                            changes))
                 #+(or)
                 (loop for slot in (distributed-slots object)
                       for slot-name = (if (symbolp slot) slot (car slot))
                       for slot-merger = (if (symbolp slot) *synchronization-merger* (cdr slot))
                       and k from 0 do
                         (unless slot-merger
                           (error "No merging defined for slot ~A of distributed object ~A"
                                  object slot-name))
                         (setf (slot-value object slot-name)
                               (funcall slot-merger
                                        object
                                        (loop for (proc . data) in changes
                                              collect (cons proc (nth k data))))))
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

(cl-store:defstore-cl-store (obj function stream)
  (cl-store::store-simple-string "<some-function>" stream))

(cl-store:defstore-cl-store (obj standard-object stream)
  (cl-store::store-simple-string (symbol-name (class-name (class-of obj))) stream))

(pushnew :ddo *features*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;  Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ddo-performance-check ()
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

(defun test-distributed ()
  (cl-mpi:mpi-init)
  (let ((test (make-number-relation 3)))
    (sb-ext:finalize test
                     (lambda ()
                       (format *trace-output* "~&Finalizing~%")
                       (force-output *trace-output*)))
    (sb-ext:gc :full t))
  (let ((test (make-number-relation 3)))
    (sb-ext:finalize test
                     (lambda ()
                       (format *trace-output* "~&Finalizing~%")
                       (force-output *trace-output*)))
    (accessing-exclusively ((table *distributed-objects*))
      (setf (gethash 1 table) test))
    (sb-ext:gc :full t)
    (accessing-exclusively ((table *distributed-objects*))
      (print (hash-table-count table)))
    ;; entry is apparently cleared only later when (gc :full t) is called separately
    )

  (let ((object (make-instance 'dummy :data 10)))
    (ensure-distributed-class (find-class 'dummy))
    (make-distributed-object object '(0))
    (sb-ext:finalize object (distributed-finalizer (local-id object)))
    (sb-ext:gc :full t)
    )
  (dbg-on :distribute)
  (reset-distributed-objects)
  (let ((object (make-instance 'dummy)))
    (make-distributed-object object '(0))
    (synchronize))
  
  (trees:pprint-tree (cdr (first (net.scipolis.relations::indices *distribution*))))
  
  (funcall (op-merger '+ 0) '((1 . 1) (2 . 2) (3 . 4)))
  )

