;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; app-utils.lisp
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

(in-package :fl.application)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun homogenized-diffusion-tensor (asv)
  "Computes the homogenized coefficient using the formula $$ A^\mathrm{hom}
= \int_Y A (\Id - \nabla N) \dy $$"
  (let* ((problem (problem asv))
	 (mesh (mesh asv))
	 (dim (dimension mesh))
	 (result (zeros dim))
	 (id-mat (eye dim))
	 (fe-class (fe-class asv)))
    (doskel (cell mesh :dimension :highest :where :surface)
      (let* ((coeffs (coefficients-of-cell cell mesh problem))
	     (diffusion-function (getf coeffs 'FL.CDR::DIFFUSION))
	     (fe (get-fe fe-class cell))
	     (qrule (quadrature-rule fe))
	     (values (transpose (get-local-from-global-vec cell fe asv))))  ; (dim x n-basis)
	(loop for shape-grads across (ip-gradients fe qrule) ; (n-basis x dim)-matrix
	      and ip in (integration-points qrule) do
	      (let* ((gradN-hat (m* values shape-grads)) ; (dim x dim)
		     (lcoords (ip-coords ip))
		     (global (local->global cell lcoords))
		     (Dphi (local->Dglobal cell lcoords))
		     (weight (* (ip-weight ip) (abs (det Dphi))))
		     (Dphi^-1 (m/ Dphi))
		     (grad-N (transpose (m* gradN-hat Dphi^-1)))
		     (diff-tensor (evaluate diffusion-function (list :global global))))
		(gemm! weight diff-tensor (m- id-mat grad-N) 1.0 result)))))
    result))

(defmethod average-coefficient (ansatz-space &key coefficient)
  (let* ((problem (problem ansatz-space))
	 (mesh (mesh ansatz-space))
	 (order (discretization-order (fe-class ansatz-space)))
	 (result nil))
    (doskel (cell mesh :dimension :highest :where :surface)
      (let* ((coeffs (coefficients-of-cell cell mesh problem))
	     (coeff-function (getf coeffs coefficient))
	     (factor-dims (mapcar #'dimension (factor-simplices cell)))
	     (qrule (gauss-rule factor-dims (1+ order))))
	(loop for ip in (integration-points qrule) do
	      (let* ((lcoords (ip-coords ip))
		     (global (local->global cell lcoords))
		     (Dphi (local->Dglobal cell lcoords))
		     (weight (* (ip-weight ip) (abs (det Dphi))))
		     (coeff-ip (evaluate coeff-function (list :global global))))
		(if result
		    (axpy! weight coeff-ip result)
		    (setq result (scal weight coeff-ip)))
		))))
    result))

(defmethod correction-tensor ((solution <ansatz-space-vector>) (rhs <ansatz-space-vector>))
  (let ((result (make-real-matrix (multiplicity solution) (multiplicity rhs))))
    (dovec ((key) solution)
      (gemm! 1.0 (vref solution key) (vref rhs key) 1.0 result :tn))
    result))

(defun convert-correction (mat)
  "Converts the (dim^2)x(dim^2) matrix returned as result of correction
tensor into an (dim x dim)-array with (dim x dim)-matrix entries."
  (let* ((dim (floor (sqrt (nrows mat))))
	 (result (make-array (list dim dim) :initial-element nil)))
    (dotimes (i dim)
      (dotimes (j dim)
	(setf (aref result i j) (make-real-matrix dim))))
    (multi-for (index (make-fixnum-vec 4) (make-fixnum-vec 4 (1- dim)))
      (setf (mref (aref result (aref index 0) (aref index 1))
		     (aref index 2) (aref index 3))
	    (mref mat
		     (+ (* dim (aref index 0)) (aref index 2))
		     (+ (* dim (aref index 1)) (aref index 3) ))))
    result))

(defun effective-tensor (blackboard)
  (with-items (&key ansatz-space problem solution rhs) blackboard
    (and solution rhs
	 (typecase problem
	   (<cdr-problem>
	    (m- (average-coefficient ansatz-space :coefficient 'FL.CDR::DIFFUSION)
		(correction-tensor solution rhs)))
	   (<elasticity-problem> 
	    (m- (average-coefficient ansatz-space :coefficient 'FL.ELASTICITY::ELASTICITY)
		(convert-correction (correction-tensor solution rhs))))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *result* nil
  "Special variable used for storing the blackbboard of the last
computation.")
