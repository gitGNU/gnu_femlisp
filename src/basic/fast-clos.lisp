;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; fast-clos.lisp
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

(in-package :macros)

;;; Better might be to define replacements of defclass and defmethod in
;;; another package...

; (in-package :cl-user)

; (defpackage "FAST-CLOS"
;   (:use "COMMON-LISP")
;   (:export ))

; (in-package :fast-clos)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; A fast interface for circumventing method dispatch, inlining methods,
;;; and making structures accessible via a defclass-like interface.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; default values
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro defgeneric* (&rest rest)
  `(defgeneric ,@rest))
(defmacro defmethod* (&rest rest)
  `(defmethod ,@rest))
(defmacro defmethod** (&rest rest)
  `(defmethod ,@rest))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; activation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro with-fast-clos (names &body body)
  "Wrap this macro around definitions to use structure-clos."
  `(macrolet
    ,(mapcar #'(lambda (name)
		 (ecase name
		   (defgeneric* `(defgeneric* (&rest rest)
				  `(defgeneric*-fast ,@rest)))
		   (defmethod* `(defmethod* (&rest rest)
				 `(defmethod*-fast ,@rest)))
		   (defmethod** `(defmethod** (&rest rest)
				  `(defmethod**-fast ,@rest)))))
	     (or names '(defgeneric* defmethod* defmethod**)))
    ,@body))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; defgeneric*
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro defgeneric*-fast (gf-name arg-list &rest rest)
  `(progn
    (c:defknown ,gf-name ,(mapcar (constantly t) arg-list) t ())
    ;; clear previous transforms
    (setf (c::function-info-transforms (ext::info function info (quote ,gf-name))) ())
    (defgeneric ,gf-name ,arg-list ,@rest)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; defmethod*
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun defmethod*-expander (gf-name arg-list body &key result-type inline)
  (let* ((arg-names (mapcar #'(lambda (arg) (if (consp arg) (car arg) arg))
			    arg-list))
	 (arg-types (mapcar #'(lambda (arg) (if (consp arg) (cadr arg) t))
			    arg-list))
	 (func-name (intern (format nil "~A ~A" gf-name arg-types))))
    (list
     'progn
     (when inline `(declaim (inline ,func-name)))
     `(declaim (ftype (function ,arg-types ,(or result-type t)) ,func-name))
     `(defun ,func-name ,arg-names ,@body)
     `(defmethod ,gf-name ,arg-list (,func-name ,@arg-names))
     `(eval-when (:load-toplevel :compile-toplevel :execute)
       (unless (c::info function c::info (quote ,gf-name))
	 (c:defknown ,gf-name ,(mapcar (constantly t) arg-names) t ())))
     `(c:deftransform ,gf-name (,arg-names ,arg-types *)
       '(,func-name ,@arg-names))
     `(ensure-generic-function ',gf-name))))

(defmacro defmethod*-fast (gf-name arg-list &body body)
  (defmethod*-expander gf-name arg-list body))

(defmacro defmethod**-fast (gf-name arg-list &body body)
  (defmethod*-expander gf-name arg-list body :inline t))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; defclass*
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmacro accelerate-slot-access (class-name slots)
  "This stuff does not yet work."
  (let ((obj (gensym "obj"))
	(val (gensym "val")))
    `(progn
      ,@(mapcan
	 #'(lambda (slot-def)
	     (destructuring-bind (slot-name &key accessor reader writer (type t) &allow-other-keys)
		 slot-def
	       (let ((slot-location
		      (pcl::slot-definition-location
		       (pcl::find-slot-definition (find-class class-name) slot-name))))
		 (when accessor
		   (transform-gf-call
		    accessor (list obj) (list class-name)
		    `(the ,type (pcl::%slot-ref (pcl::std-instance-slots obj) ,slot-location))))
		 (when accessor
		   (transform-gf-call
		    `(setf ,accessor) (list val obj) (list type class-name)
		    `(setf (pcl::%slot-ref (pcl::std-instance-slots obj) ,slot-location) ,val)))
		 (when reader
		   (transform-gf-call
		    reader (list obj) (list class-name)
		    `(the ,type (pcl::%slot-ref (pcl::std-instance-slots obj) ,slot-location))))
		 (when writer
		   (transform-gf-call
		    writer (list val obj) (list type class-name)
		    `(setf (pcl::%slot-ref (pcl::std-instance-slots obj) ,slot-location) ,val))))))
	 slots))))

;;;; Testing: (test-structure-clos)

(with-fast-clos ()
  (defmethod** my-elt ((vec array) i) (aref vec i))
  (defmethod** my-elt ((list list) i) (nth i list)))

(defun test-structure-clos ()
  (princ (c::function-info-or-lose 'my-elt))
  (time
   (let ((vec (make-array 10 :element-type 'double-float :initial-element 0.0d0)))
     (declare (optimize (speed 3) (safety 1)))
     (declare (type (simple-array double-float (*)) vec))
     (let ((sum 0.0d0))
       (declare (type double-float sum))
       (dotimes (k 100000 sum)
	 (dotimes (i 10)
	   (incf sum (my-elt vec i)))))))
)


