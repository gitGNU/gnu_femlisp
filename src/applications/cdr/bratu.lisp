;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; bratu.lisp - Solve the Bratu problem
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

(in-package :application)

(defvar *u_1/2-observe*
  (list
   (list (format nil "~19@A" "u(midpoint)") "~19,10,2E"
	 #'(lambda (blackboard)
	     (let* ((sol (getbb blackboard :solution))
		    (dim (dimension (mesh (ansatz-space sol))))
		    (val (fe-value sol (make-double-vec dim 0.5))))
	       (vref (aref val 0) 0))))))

(defun bratu-computation (dim &key (plot t) (time 5.0) (output 1))
  "bratu-~Dd - Solve the Bratu problem in ~DD

Solves the Bratu problem -Delta u +e^u =0.  It reports the value in the
midpoint of the domain."
  (defparameter *result*
    (solve (blackboard
	    :problem (bratu-problem dim)
	    :success-if `(> :time ,time) :output output
	    :observe (append *stationary-fe-strategy-observe* *u_1/2-observe*))))
  (when plot (plot (getbb *result* :solution))))

(defun make-bratu-demo (dim)
  (multiple-value-bind (title short long)
      (extract-demo-strings (documentation 'bratu-computation 'function))
    (let ((demo
	   (make-demo
	    :name (format nil title dim)
	    :short (format nil short dim)
	    :long long
	    :execute (lambda () (bratu-computation dim :plot t :output 1)))))
      (adjoin-demo demo *cdr-demo*))))

;;; 2D and 3D Bratu problem
(make-bratu-demo 1)
(make-bratu-demo 2)
(make-bratu-demo 3)

(defun bratu-tests ()
  (bratu-computation 1))

;;; (application::bratu-tests)
(fl.tests::adjoin-test 'bratu-tests)
