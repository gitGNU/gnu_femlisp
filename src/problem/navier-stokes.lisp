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
(defpackage "FL.NAVIER-STOKES"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.ALGEBRA" "FL.FUNCTION" "FL.MESH" "FL.PROBLEM")
  (:export
   "<NAVIER-STOKES-PROBLEM>" "NO-SLIP-BOUNDARY"
   "STANDARD-NAVIER-STOKES-PROBLEM" "UNIT-VECTOR-FORCE"
   "DRIVEN-CAVITY" "PERIODIC-CAVITY")
  (:documentation "Defines the class of Navier-Stokes problems."))
(in-package "FL.NAVIER-STOKES")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <navier-stokes-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <navier-stokes-problem> (<pde-problem>)
  ()
  (:documentation "Navier-Stokes problem."))

(defmethod nr-of-components ((problem <navier-stokes-problem>))
  (1+ (dimension (domain problem))))

(defmethod interior-coefficients ((problem <navier-stokes-problem>))
  "Interior coefficients for the Navier-Stokes problem."
  '(VISCOSITY REYNOLDS FORCE))

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
	(values (make-array (1+ dim) :initial-element #m(0.0))))
    (setf (aref flags dim) nil)
  (constant-coefficient flags values)))

(defun unit-vector-force (dim &optional (direction 0))
  (constant-coefficient
   (let ((result (make-array (1+ dim) :initial-element #m((0.0)))))
     (setf (aref result direction) #m((1.0)))
     result)))

(defun standard-navier-stokes-problem (domain &key (viscosity 1.0) (reynolds 0.0) force)
  (let ((dim (dimension domain)))
    (make-instance
     '<navier-stokes-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (cond ((member-of-skeleton? patch (domain-boundary domain))
		(list 'FL.NAVIER-STOKES::CONSTRAINT (no-slip-boundary dim)))
	       ((= dim (dimension patch))
		(list 'FL.NAVIER-STOKES::VISCOSITY (ensure-coefficient viscosity)
		      'FL.NAVIER-STOKES::REYNOLDS (ensure-coefficient reynolds)
		      'FL.NAVIER-STOKES::FORCE (ensure-coefficient force)))
	       (t ()))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Driven cavity
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun driven-cavity-upper-boundary (dim &optional smooth-p)
  (let ((flags (make-array (1+ dim) :initial-element (not smooth-p)))
	(values (make-array (1+ dim) :initial-element #m(0.0))))
    (setf (aref flags dim) nil) ; no constraint for pressure
    (when smooth-p  ; zero constraint for u_n
      (setf (aref flags (1- dim)) t))
    (unless smooth-p
      (setf (aref values 0) #m(1.0)))
    (constant-coefficient flags values)))

(defun driven-cavity-force (dim)
  (unit-vector-force dim 0))

(defun driven-cavity (dim &key (viscosity 1.0) (reynolds 0.0) smooth-p)
  (let ((domain (n-cube-domain dim)))
    (make-instance
     '<navier-stokes-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (let ((midpoint (midpoint patch)))
	   (cond
	     ((= dim (dimension patch)) ; interior coefficients
	      (list 'FL.NAVIER-STOKES::VISCOSITY (ensure-coefficient viscosity)
		    'FL.NAVIER-STOKES::REYNOLDS (ensure-coefficient reynolds)))
	      ;; upper boundary
	     ((and (= (dimension patch) (1- dim))
		   (= (aref midpoint (1- dim)) 1.0))
	      (append (when smooth-p
			(list 'FL.NAVIER-STOKES::FORCE (driven-cavity-force dim)))
		      (list 'FL.NAVIER-STOKES::CONSTRAINT
			    (driven-cavity-upper-boundary dim smooth-p))))
	     ;; lower corner sets also pressure to zero:
	     ((mzerop midpoint)
	      (list 'FL.NAVIER-STOKES::CONSTRAINT
		    (constraint-coefficient (1+ dim) 1)))
	     ;; other boundaries have standard no-slip bc
	     (t (list 'FL.NAVIER-STOKES::CONSTRAINT
		      (no-slip-boundary dim))))))
     :properties (list :linear-p (zerop reynolds)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Periodic cavity
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun oscillating-force (dim)
  (function->coefficient
   #'(lambda (x)
    (let ((result (make-array (1+ dim))))
      (dotimes (i (1+ dim) result)
	(setf (aref result i)
	      (if (< i dim)
		  (make-real-matrix `((,#I"sin(2*pi*x[0])*cos(2*pi*x[1])")))
		  (zeros 1))))))))

(defun periodic-cavity (dim &key (viscosity 1.0) (reynolds 0.0))
  (let ((domain (n-cell-domain dim)))
    (make-instance
     '<navier-stokes-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (when (= (dimension patch) dim)
	   (list 'VISCOSITY (ensure-coefficient viscosity)
		 'REYNOLDS (ensure-coefficient reynolds)
		 'FORCE (oscillating-force dim)
		 ))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Channel flow
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun simple-pressure-boundary-conditions (dim dir value)
  "This function is a preliminary version which works only for boundaries
lying in a coordinate hyperplane."
  (assert (<= 0 dir (1- dim)))
  (let ((flags (make-array (1+ dim) :initial-element t))
	(values (make-array (1+ dim) :initial-element #m(0.0))))
    (setf (aref flags dir) nil) ; no constraint for dir-component
    (setf (aref values dim) (ensure-matlisp value))
    (constant-coefficient flags values)))

(defun cubic-channel (dim &key (viscosity 1.0) (reynolds 0.0) (direction 0))
  "The channel is a simple test problem which has for every Reynolds number
a Hagen-Poiseuille solution with linear pressure.  The corresponding
velocity profile is a solution to @math{-\Delta u = constant}."
  (let ((domain (n-cube-domain dim)))
    (make-instance
     '<navier-stokes-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (let ((dir-coordinate (aref (midpoint patch) direction)))
	   (cond
	     ((= dim (dimension patch)) ; interior coefficients
	      (list 'FL.NAVIER-STOKES::VISCOSITY (ensure-coefficient viscosity)
		    'FL.NAVIER-STOKES::REYNOLDS (ensure-coefficient reynolds)))
	     ((and (= (1- dim) (dimension patch))
		   (not (/= dir-coordinate 0.0 1.0)))
	      (list 'FL.NAVIER-STOKES::CONSTRAINT
		    (simple-pressure-boundary-conditions
		     dim direction (if (zerop dir-coordinate) 1.0 0.0))))
	     ;; other boundaries have standard no-slip bc
	     (t (list 'FL.NAVIER-STOKES::CONSTRAINT (no-slip-boundary dim))))))
     :properties (list :linear-p (zerop reynolds)))))

;;; Testing: (test-navier-stokes)

(defun test-navier-stokes ()
  (describe
   (standard-navier-stokes-problem
    (n-cube-domain 2) :force (constantly (unit-vector 2 0))))
  (describe (driven-cavity 2 :smooth-p nil :reynolds 1.0))
  )

(fl.tests:adjoin-test 'test-navier-stokes)


