;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; newton.lisp
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

;;;; Purpose: nonlinear solver

;;;; This is old code.  It was used only for odes in a NM1 course.

(in-package :iteration)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; newton iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun newton (f x-start &key output
	       approximate-gradient approximate-gradient-inverse
	       (approximate-gradient-inverter (constantly #'identity))
	       (threshold 1.0d-12) (maxsteps 100))
  "This function solves the equation f(x)=0 by an approximate Newton's
method. x-start is the initial value.  Usually, you should also
provide EITHER an approximation of the gradient of f as a function
giving a number or matrix OR an approximation to the inverse of this
gradient OR an approximate inverse of the gradient as a mapping
between vectors.  If you don't, then the method defaults to a
fixed-point iteration for $identity-f$."
  (when output (format t "Newton iteration:~%"))
  (loop for x = x-start
	then (m- x (cond (approximate-gradient
			  (m/ (funcall approximate-gradient x) f-value))
			 (approximate-gradient-inverse
			  (m* (funcall approximate-gradient-inverse x) f-value))
			 (approximate-gradient-inverter
			  (funcall (funcall approximate-gradient-inverter x)
				   f-value))))
	for i = 0 then (1+ i) until (> i maxsteps)
	for f-value = (funcall f x)
	when output do (format t "Step ~3D: |f| = ~8,3,2E, x = ~A~%" i (norm f-value) x)
	until (< (norm f-value) threshold)
	finally (return x)))

(defun test-newton ()
  (flet ((simple-linear-f (x) (+ 1.0 (* 0.5 x))))
    (newton #'simple-linear-f 1.0 :output t)
    (newton #'simple-linear-f 1.0 :approximate-gradient (constantly 0.5) :output t))
  )

  
;;;; Testing: (test-newton)
(adjoin-femlisp-test 'test-newton)