;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; elasticity.lisp
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

(defpackage "FL.ELASTICITY"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"FL.ALGEBRA" "FL.FUNCTION" "FL.MESH" "FL.PROBLEM")
  (:export
   "<ELASTICITY-PROBLEM>" "ISOTROPIC-ELASTICITY-TENSOR"
   "CHECK-ELASTICITY-TENSOR"
   "SYSTEM-DIFFUSION-PROBLEM" "STANDARD-ELASTICITY-PROBLEM")
  (:documentation "Defines elasticity problems."))

(in-package :fl.elasticity)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <elasticity-problem>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <elasticity-problem> (<pde-problem>)
  ()
  (:documentation "Elasticity problems."))

(defmethod nr-of-components ((problem <elasticity-problem>))
  (dimension (domain problem)))

(defmethod interior-coefficients ((problem <elasticity-problem>))
  "Interior coefficients for the elasticity problem."
  '(ELASTICITY FORCE GAMMA))
(defmethod boundary-coefficients ((problem <elasticity-problem>))
  "Boundary coefficients for the elasticity problem."
  '(CONSTRAINT))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of standard elasticity problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun isotropic-elasticity-tensor (&key dim lambda mu)
  "Returns the tensor corresponding to the Lam'e constants @math{lambda}
and @math{mu}, i.e.: @math{A_{ij}^{kl} = lambda delta_{ik} delta_{jl} + mu
(delta_{ij} delta_{kl} + delta_{kj} delta_{il})}."
  (let ((tensor (make-array (list dim dim) :initial-element nil)))
    (dotimes (k dim)
      (dotimes (l dim)
	(let ((mat (make-real-matrix dim)))
	  (dotimes (i dim)
	    (dotimes (j dim)
	      (setf (mref mat i j)
		    (+ (* lambda (if (and (= i k) (= j l)) 1.0 0.0))
		       (* mu  (+ (if (and (= i j) (= k l)) 1.0 0.0)
				 (if (and (= k j) (= i l)) 1.0 0.0)))))))
	  (setf (aref tensor k l) mat))))
    tensor))

(defun check-elasticity-tensor (tensor &optional (threshold 1.0e-6))
  "Checks the symmetries in the elasticity tensor."
  (labels ((tref (index)
	     (mref (aref tensor (aref index 0) (aref index 1))
		   (aref index 2) (aref index 3)))
	   (same? (ind1 ind2)
	     (<= (abs (- (tref ind1) (tref ind2))) threshold)))
    (let ((dim (array-dimension tensor 0)))
      (multi-for (index (make-fixnum-vec 4) (make-fixnum-vec 4 (1- dim)))
	(flet ((check-permutation (permutation)
		 (let ((permuted (permute permutation index)))
		   (unless (same? index permuted)
		     (format t "~A: t[~A] = ~A but t[~A] = ~A~%" permutation
			     index (tref index) permuted (tref permuted))))))
	  (check-permutation #(1 0 3 2))
	  (check-permutation #(2 1 0 3)))))
    tensor))

(defun standard-elasticity-problem (domain &key lambda mu force)
  (let ((dim (dimension domain)))
    (make-instance
     '<elasticity-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (if (member-of-skeleton? patch (domain-boundary domain))
	     (list 'CONSTRAINT (constraint-coefficient dim 1))
	     (list 'ELASTICITY
		   (constant-coefficient
		    (isotropic-elasticity-tensor :dim dim :lambda lambda :mu mu))
		   'FORCE
		   (ensure-coefficient force)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Systems of diffusion equations (mainly for testing purposes)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun system-diffusion-tensor (&key dim nr-comps D)
  "Returns a tensor for a system of separate diffusion equations."
  (let ((tensor (make-array (list nr-comps nr-comps) :initial-element nil)))
    (dotimes (k nr-comps)
      (dotimes (l nr-comps)
	(setf (aref tensor k l)
	      (if (= k l)
		  (scal! D (eye dim))
		  (zeros dim)))))
    tensor))

(defun system-diffusion-problem (domain &key nr-comps D force)
  (let ((dim (dimension domain)))
    (make-instance
     '<elasticity-problem>
     :domain domain
     :patch->coefficients
     #'(lambda (patch)
	 (if (member-of-skeleton? patch (domain-boundary domain))
	     (list 'CONSTRAINT (constraint-coefficient nr-comps 1))
	     (list 'ELASTICITY
		   (constant-coefficient
		    (system-diffusion-tensor :dim dim :nr-comps nr-comps :D D))
		   'FORCE
		   (ensure-coefficient force)))))))


;;; Testing: (test-elasticity)

(defun test-elasticity ()
  (check-elasticity-tensor
   (isotropic-elasticity-tensor :dim 2 :lambda 1.0 :mu 2.0))
  ;; the following should be equal to (cdr-model-problem 1)
  (let ((dim 1))
    (standard-elasticity-problem
     (n-cube-domain dim) :lambda 1.0 :mu 1.0
     :force (constant-coefficient (make-array dim :initial-element (zeros 1)))))
  )

(fl.tests:adjoin-test 'test-elasticity)


