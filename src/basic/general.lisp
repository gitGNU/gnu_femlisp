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

(defun runtime-compile (source)
  "Calls compile on the provided @arg{source}.  When :compile is activated
for debugging, the source code is printed."
  (let ((*print-circle* nil))
    (dbg :compile "Compiling source:~%~A" source))
  (compile nil source))

