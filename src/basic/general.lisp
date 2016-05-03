;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; general.lisp - globally useful definitions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2006- Nicolas Neuss, University of Karlsruhe.
;;; All rights reserved.
;;; 
;;; Redistribution and use in source and binary forms, with or without
;;; modification, are permitted provided that the following conditions are
;;; met:
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
;;; MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITIES HEIDELBERG AND KARLSRUHE OR
;;; OTHER CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;; SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;;; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;;; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;;; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
;;; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.utilities)

;;;; This files provides definitions which have influence across several
;;;; other packages.

(defgeneric make-analog (obj)
  (:documentation "Generate an analogous but empty data structure.")
  (:method (obj) (make-instance (class-of obj))))

(defun file-documentation (docstring)
  "If the manual is sorted by file, the string handed to this function
describes the use of the respective file."
  (declare (ignore docstring)))

(defun concept-documentation (docstring)
  "Documents a certain concept."
  (declare (ignore docstring)))

(defclass property-mixin ()
  ((properties :accessor properties :initform () :initarg :properties
	       :type list :documentation "A property list for storing
unstructured information about this object."))
  (:documentation "A mixin adding a property slot to the class."))

;;; old interface

(defun get-property (object property)
  "Gets @arg{property} for @arg{object}.  Returns NIL also if
@arg{property} is not available."
  (getf (slot-value object 'properties) property))

(defun (setf get-property) (value object property)
  "Sets the property @arg{property} of @arg{problem} to @arg{value}."
  (setf (getf (slot-value object 'properties) property) value))

;;; new interface using CLOS (WITH-SLOTS/SLOT-VALUE work)
(defmethod shared-initialize :after ((object property-mixin) slot-names
				     &rest initargs &key &allow-other-keys)
  (declare (ignore slot-names))
  (let* ((class (class-of object))
	 (ordinary-slot-initargs (mappend #'fl.amop:slot-definition-initargs
					  (fl.amop:compute-slots class))))
    (loop for (key value) on initargs by #'cddr
	 unless (member key ordinary-slot-initargs) do
	 (setf (slot-value object (intern (symbol-name key) *package*))
	       value))))

(defmethod slot-missing (class (object property-mixin) slot-name
			 (operation (eql 'slot-value)) &optional new-value)
  (declare (ignore class new-value))
  (getf (slot-value object 'properties) slot-name))

(defmethod slot-missing (class (object property-mixin) slot-name
			 (operation (eql 'setf)) &optional new-value)
  (declare (ignore class))
  (setf (getf (slot-value object 'properties) slot-name) new-value))

(defmethod slot-missing (class (object property-mixin) slot-name
			 (operation (eql 'slot-boundp)) &optional new-value)
  (declare (ignore class slot-name new-value))
  t)

(defmethod slot-missing (class (object property-mixin) slot-name
			 (operation (eql 'slot-makunbound)) &optional new-value)
  (declare (ignore class new-value))
  (error "Property slot '~A' cannot be unbound" slot-name))


;;;; Result modification by hooks

(defparameter *hooks*
  (make-hash-table)
  "Table mapping a symbol to a list of hooks.
A hook is a function which may be called on a property object
for updating other properties.")

(defun call-hooks (function-name object)
  "Call all hooks defined for @arg{function-name} on @arg{object} and
returns @arg{object}."
  (loop for hook-entry in (gethash function-name *hooks*) do
    (funcall (cdr hook-entry) object))
  object)

(defun add-hook (function-name hook-name hook-function)
  "Add a hook with name @arg{hook-name} to the hooks for @arg{function-name}."
  (let ((hook-list (gethash function-name *hooks*)))
    (if (null hook-function)
        (setf hook-list (remove hook-name hook-list :key #'car))
        (let ((entry (assoc hook-name hook-list)))
          (if entry
              ;; replace existing setting
              (setf (cdr entry) hook-function)
              (pushnew (cons hook-name hook-function) hook-list))))
    (setf (gethash function-name *hooks*) hook-list)))

;;;; Testing

(defun test-general ()
  (loop+ ((i (range :below 0))) do (error "should not be reached: ~D" i))
  (let ((x (make-array 10 :initial-element 0))
	(y (make-list 10 :initial-element 1)))
    (loop+ ((xc x) (yc y) i) doing
       (setf xc (+ i yc))
       finally (return x)))
  (loop+ ((k '(1 2 3)) i) collecting (+ k i))

  #-cmu
  (let ((x (make-instance 'property-mixin :properties (list :a 1))))
    (assert (= (slot-value x :a) 1))
    (setf (slot-value x :a) 3)
    (assert (= (slot-value x :a) 3)))

  (let ((x (make-instance 'property-mixin :properties '(:a 1)))
	(n 100))
    (time (loop repeat n do (get-property x :a)))
    #-cmu (time (loop repeat n do (slot-value x :a)))
    )
  (let ((x (make-instance 'property-mixin)))
    (describe x)
    (with-slots (x) x
      x))
  )

;;; (test-general)
(fl.tests:adjoin-test 'test-general)
