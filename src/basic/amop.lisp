;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; amop.lisp - AMOP extensions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
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

(defpackage "FL.AMOP"
  (:use "COMMON-LISP"
	#+CMU "MOP" #+SBCL "SB-MOP")
  (:export "FIND-PROGRAMMATIC-CLASS" "MAKE-PROGRAMMATIC-INSTANCE"
	   "REMOVE-SUBCLASS-METHODS" "CLASS-DIRECT-SUPERCLASSES")
  (:documentation
   "This package provides some MOP functionality.  These functions are
non-ANSI and may represent a problem when porting Femlisp."))

(in-package :fl.amop)

(defun find-programmatic-class (superclasses)
  "Finds and, if necessary, generates a class from the given superclasses."
  (cond
    ((null superclasses) T)
    ((null (cdr superclasses)) (car superclasses))
    (t (let ((class (find-if
		     #'(lambda (class)
			 (equal superclasses
				(class-direct-superclasses class)))
		     (class-direct-subclasses (car superclasses)))))
	 (or class
	     (make-instance 'standard-class
			    :name (intern (format nil "~A" (mapcar #'class-name superclasses)))
			    :direct-superclasses superclasses
			    :direct-slots ()))))))

(defun make-programmatic-instance (superclass-names &rest initargs)
  "Makes an instance of a class denoted by a list of the names of its
superclasses.  This class is generated automatically, if necessary."
  (apply #'make-instance
	 (if (symbolp superclass-names)
	     (find-class superclass-names)
	     (find-programmatic-class (mapcar #'find-class superclass-names)))
         initargs))

(defun remove-subclass-methods (gf template-args)
  "Removes all methods dispatching on subclasses of the template
arguments."
  (loop for method in (copy-seq (generic-function-methods gf))
	when (every #'subtypep
		    (method-specializers method)
		    (mapcar #'(lambda (arg)
				(if (consp arg)
				    (second arg)
				    T))
			    template-args))
	do (remove-method gf method)))
