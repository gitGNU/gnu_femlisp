;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; general.lisp - globally useful definitions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003 Nicolas Neuss, University of Heidelberg.
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
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF HEIDELBERG OR OTHER
;;; CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
;;; EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
;;; PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
;;; PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
;;; LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
;;; NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
;;; SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(in-package :fl.utilities)

;;;; This files provides definitions which have influence across several
;;;; other packages.

(defgeneric make-analog (obj)
  (:documentation "Generate an analogous but empty data structure."))

(defun file-documentation (docstring)
  "If the manual is sorted by file, the string handed to this function
describes the use of the respective file."
  (declare (ignore docstring)))

(defclass property-mixin ()
  ((properties :accessor properties :initform () :initarg :properties
	       :type list :documentation
    "A property list which is used to store unstructured information about
this object."))
  (:documentation "A mixin which adds a slot of properties to the class."))

(defun property-set-p (object property)
  "Returns T if @arg{property} is found in the object's properties."
  (get-properties (slot-value object 'properties) (list property)))

(defun get-property (object property)
  "Gets @arg{property} for @arg{object}."
  (getf (slot-value object 'properties) property))

(defun (setf get-property) (value object property)
  "Sets the property @arg{property} of @arg{problem} to @arg{value}."
  (setf (getf (slot-value object 'properties) property) value))

(defmacro with-properties (properties object &body body)
  "Work with @arg{properties} on the property list of @arg{object}."
  (with-gensyms (obj)
    `(let ((,obj ,object))
       (symbol-macrolet ,(loop for prop in properties collect
			      `(,prop (getf (properties ,obj) ',prop)))
	 ,@body))))


(defun runtime-compile (source)
  "Calls compile on the provided @arg{source}.  When :compile is activated
for debugging, the source code is printed."
  (let ((*print-circle* nil))
    (dbg :compile "Compiling source:~%~A" source))
  (compile nil source))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass range ()
  ((from :initform 0 :initarg :from)
   (to :initform nil :initarg :to)
   (below :initform nil :initarg :below)
   (by :initform 1 :initarg :by))
  (:documentation "Range of numbers for iteration."))

(defun range (&rest args &key to below)
  "Constructor for a range of numbers."
  (when (and to below)
    (error "Only one limit should be given."))
  (apply #'make-instance 'range args))

(defgeneric iterator (vec)
  (:documentation "Returns an iterator for the given vector.")
  (:method ((range range)) (slot-value range 'from))
  (:method ((vec vector)) 0)
  (:method ((vec list)) vec))

(defgeneric iterator-next (vec iterator)
  (:documentation "Returns an incremented @arg{iterator}.")
  (:method ((range range) i) (+ i (slot-value range 'by)))
  (:method ((vec vector) i) (1+ i))
  (:method ((vec list) tail) (cdr tail)))

(defgeneric iterator-end-p (vec iterator)
  (:method ((range range) i)
    (with-slots (from to below) range
      (cond (to (> i to))
	    (below (>= i below)))))
  (:method ((vec vector) i) (>= i (length vec)))
  (:method ((vec list) tail) (null tail)))
  
(defgeneric reference (vec iterator)
  (:documentation "Reader for the element of @arg{vec} referenced by @arg{iterator}.")
  (:method ((range range) i) i)
  (:method ((vec vector) i) (aref vec i))
  (:method ((vec list) tail) (car tail)))
  
(defgeneric (setf reference) (value vec iterator)
  (:documentation "Setter for the element of @arg{vec} referenced by @arg{iterator}.")
  (:method (value (vec vector) i)
    (setf (aref vec i) value))
  (:method (value (vec list) tail)
    (setf (car tail) value)))

(defmacro loop+ (items &body body)
  "Iterates @arg{body} over @arg{items}."
  (let ((vectors (loop for i below (length items)
		      collect (gensym (format nil "V~D" i))))
	(iterators (loop for i below (length items)
		      collect (gensym (format nil "I~D" i)))))
    `(let ,(mapcar #'(lambda (vector item)
		       (list vector (if (listp item) (cadr item) `(range))))
		   vectors items)
       (symbol-macrolet ,(mapcar #'(lambda (vector iterator item)
				     `(,(if (listp item) (car item) item)
					(reference ,vector ,iterator)))
				 vectors iterators items)
	   (loop ,@(loop for iterator in iterators and vector in vectors appending
			 `(for ,iterator = (iterator ,vector)
			   then (iterator-next ,vector ,iterator)))
		 ,@(loop for iterator in iterators and vector in vectors appending
			 `(until (iterator-end-p ,vector ,iterator)))
	    ,@body
	 )))))

;;;; Testing

(defun test-general ()
  (loop+ ((i (range :below 0))) do (error "should not be reached"))
  (let ((x (make-array 10))
	(y (make-list 10 :initial-element 1)))
    (loop+ ((xc x) (yc y) i) doing
       (setf xc (+ i yc))
       finally (return x)))
  )

;;; (test-general)
(fl.tests:adjoin-test 'test-general)