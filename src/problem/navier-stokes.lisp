;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; navier-stokes.lisp
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

(in-package :cl-user)
(defpackage "NAVIER-STOKES"
  (:use "COMMON-LISP" "MATLISP" "MACROS" "UTILITIES" "ALGEBRA" "MESH" "PROBLEM")
  (:export
   "<NAVIER-STOKES-PROBLEM>" "NO-SLIP-BOUNDARY"
   "STANDARD-NAVIER-STOKES-PROBLEM" "UNIT-VECTOR-FORCE"
   "DRIVEN-CAVITY" "PERIODIC-CAVITY"))
(in-package :navier-stokes)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <navier-stokes-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <navier-stokes-problem> (<problem>)
  ()
  (:documentation "Navier-Stokes problem."))

(defmethod interior-coefficients ((problem <navier-stokes-problem>))
  "Interior coefficients for the Navier-Stokes problem."
  '(VISCOSITY REYNOLDS FORCE))
(defmethod boundary-coefficients ((problem <navier-stokes-problem>))
  "Boundary coefficients for the Navier-Stokes problem."
  '(CONSTRAINT))

;;; Warning: in elasticity-fe we handle incompressible Navier-Stokes with a
;;; viscosity term of the form nu*Delta u.  For non-Dirichlet boundary
;;; conditions for the velocity one should use instead a term of the form
;;; div(nu*eps(u)) with eps denoting the strain tensor.  We will switch to
;;; this alternative form soon.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of a standard Navier-Stokes problem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun no-slip-boundary (dim)
  (let ((flags (make-array (1+ dim) :initial-element t))
	(values (make-double-vec (1+ dim))))
    (setf (aref flags dim) nil)
  (constant-coefficient flags values)))

(defun unit-vector-force (dim &optional (direction 0))
  (constant-coefficient
   (let ((result (make-array (1+ dim) :initial-element [0.0])))
     (setf (aref result direction) [1.0])
     result)))

(defun standard-navier-stokes-problem (domain &key (viscosity 1.0) (reynolds 0.0) force)
  (let ((dim (dimension domain)))
    (make-instance
     '<navier-stokes-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (cond ((member-of-skeleton? patch (domain-boundary domain))
		(list 'CONSTRAINT (no-slip-boundary dim)))
	       ((= dim (dimension patch))
		(list 'VISCOSITY (ensure-coefficient viscosity)
		      'REYNOLDS (ensure-coefficient reynolds)
		      'FORCE (ensure-coefficient force)))
	       (t ()))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Driven cavity
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun driven-cavity-upper-boundary (dim)
  (let ((flags (make-array (1+ dim) :initial-element t))
	(values (make-double-vec (1+ dim))))
    (setf (aref flags dim) nil) ; no constraint for pressure
    (setf (aref values 0) 1.0d0)
    (constant-coefficient flags values)))

(defun driven-cavity-force (dim)
  (unit-vector-force dim 0))

(defun driven-cavity (dim &key (viscosity 1.0) (reynolds 0.0) smooth-p)
  (let* ((domain (n-cube-domain dim))
	 (upper (find-cell
		 #'(lambda (cell)
		     (and (= (dimension cell) (1- dim))
			  (= (aref (midpoint cell) (1- dim)) 1.0)))
		 domain)))
    (assert upper)
    (make-instance
     '<navier-stokes-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (cond
	   ((member-of-skeleton? patch (domain-boundary domain))
	    ;; boundary coeffs
	    (if (eq patch upper)
		(if smooth-p
		    (list 'FORCE (driven-cavity-force dim))
		    (list 'CONSTRAINT (driven-cavity-upper-boundary dim)))
		(list 'CONSTRAINT (no-slip-boundary dim))))
	   ;; inner coeffs
	   ((= dim (dimension patch))
	    (list 'VISCOSITY (constant-coefficient viscosity)
		  'REYNOLDS (constant-coefficient reynolds)))
	   (t ()))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Periodic cavity
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun periodic-cavity-force (dim)
  (function->coefficient
   #'(lambda (x)
    (coerce (append (make-list dim :initial-element [#I"sin(2*pi*x[0])*cos(2*pi*x[1])"])
		    (list [0.0]))
	    'vector))))

(defun periodic-cavity (dim &key (viscosity 1.0) (reynolds 0.0))
  (let ((domain (n-cell-domain dim)))
    (make-instance
     '<navier-stokes-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (when (= (dimension patch) dim)
	   (list 'VISCOSITY (constant-coefficient viscosity)
		 'REYNOLDS (constant-coefficient reynolds)
		 'FORCE (periodic-cavity-force dim)
		 ))))))



;;; Testing: (test-navier-stokes)

(defun test-navier-stokes ()
  (problem-info
   (standard-navier-stokes-problem
    *unit-quadrangle-domain* :force (constantly (unit-vector 2 0))))
  (problem-info (driven-cavity 2 :smooth-p nil))
  )

(tests:adjoin-femlisp-test 'test-navier-stokes)


