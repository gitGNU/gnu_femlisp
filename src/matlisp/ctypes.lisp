;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ctypes.lisp - Some C data types and uniform arrays
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003, 2004 Nicolas Neuss, University of Heidelberg.
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


(in-package :fl.matlisp)

(file-documentation
 "Introduce some data types for interfacing with foreign code.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; double-vec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype double-vec ()
  "Uniform @type{double-float} vector."
  '(simple-array double-float (*)))

(definline double-vec (&rest comps)
  "Returns a @class{double-vec} with the entries in @arg{comps}."
  (coerce comps 'double-vec))

(definline make-double-vec (dim &optional (init 0.0))
  "Returns a @class{double-vec} of length @arg{dim} and initial value
@arg{init}."
  (make-array dim :element-type 'double-float
	      :initial-element (float init 0.0)))

(eval-when (:compile-toplevel :load-toplevel :execute)
  (set-dispatch-macro-character
   #\# #\d  ; dispatch on #d for double-vec
   #'(lambda (stream char n)
       (declare (ignore char n))
       (let ((list (read stream nil (values) t)))
	 `(coerce ',list 'double-vec)))))

(defun unit-vector (dim i)
  "Returns a freshly created copy of the @arg{i}-th carthesian unit vector
in dimension @arg{dim}."
  (let ((vec (make-double-vec dim)))
    (setf (aref vec i) 1.0)
    vec))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; int-vec, uint-vec
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype int () '(signed-byte 32))
(deftype uint () '(unsigned-byte 32))

(deftype int-vec ()
  "Uniform @type{int} vector."
  '(simple-array int (*)))
(deftype uint-vec ()
  "Uniform @type{uint} vector."
  '(simple-array uint (*)))

(defun make-int-vec (dim &optional (init 0))
  (make-array dim :element-type 'int :initial-element init))
(defun make-uint-vec (dim &optional (init 0))
  (make-array dim :element-type 'uint :initial-element init))
