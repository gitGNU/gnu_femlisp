;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  navier-stokes-fe.lisp
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
(defpackage "NAVIER-STOKES-FE"
  (:use "COMMON-LISP" "MATLISP" "MACROS" "UTILITIES" "ALGEBRA"
	"MESH" "PROBLEM" "NAVIER-STOKES" "DISCRETIZATION" "NAVIER-STOKES")
  (:export "NAVIER-STOKES-LAGRANGE-FE" "BOUNDARY-LAST-ORDER"))
(in-package :navier-stokes-fe)

(defclass <navier-stokes-lagrange-fe> (<vector-fe-discretization>)
  ()
  (:documentation "Up to now (Q_{k+1})^dim/Q_k finite elements."))

(defun navier-stokes-lagrange-fe (order dim delta)
  "Returns (Q_{k+delta})^dim/Q_k finite elements."
  (make-instance
   '<navier-stokes-lagrange-fe>
   :components
   (append
    (make-list dim :initial-element (lagrange-fe (+ order delta)))
    (list (lagrange-fe order :integration-order (+ order delta))))))

(defmethod boundary-last-order ((asa <ansatz-space-automorphism>))
  "To solve the arising singular system with a direct decomposition we
reorder the equations such that the problem appears only at the end.  Then
rounding errors make it solvable.  Of course, an iterative solution should
be preferred."
  (let* ((mesh (mesh asa))
	 (boundary (domain-boundary (domain mesh)))
	 (bdry-keys ())
	 (inner-keys ()))
    (for-each-row-key
     #'(lambda (key)
	 (let ((cell (representative key)))
	   (if (member-of-skeleton? (patch-of-cell cell mesh) boundary)
	       (push cell bdry-keys)
	       (push cell inner-keys))))
     asa)
    (append inner-keys bdry-keys)))

(defun pressure-diagonal-inverter (mat)
  "For use in a direct decomposition or ILU smoother.  Does not really work
up to now."
  (let ((copy-mat (copy mat))
	(dim-1 (1- (nrows mat))))
    (multiple-value-bind (lr pivot info)
	(getrf! copy-mat)
      (when (or (not (eq info t))
		(> (norm copy-mat) (* 1e10 (algebra::det-from-lr lr pivot))))
	(copy! mat copy-mat)
	(setf (mat-ref copy-mat dim-1 dim-1) 1.0d0)
	(multiple-value-setq (lr pivot) (getrf! copy-mat)))
      (getrs! lr pivot (eye (nrows mat))))))


(defmethod discretize-locally ((problem <navier-stokes-problem>)
			       coeffs vecfe qrule fe-geometry
			       &key local-mat local-rhs local-sol local-u local-v)
  "Local discretization for a Navier-Stokes problem."
  (declare (type (array real-matrix (* *)) local-mat)
	   (type (or null (array real-matrix (*)))
		 local-sol local-rhs local-u local-v))
  ;; and this situation should be checked before use
  (assert (and (null local-u) (null local-v)))
  
  (let* ((nr-comps (nr-of-components vecfe))
	 (dim (1- nr-comps))
	 (cell-dim (dimension (getf fe-geometry :cell)))
	 (viscosity (getf coeffs 'NAVIER-STOKES::VISCOSITY))
	 (reynolds (getf coeffs 'NAVIER-STOKES::REYNOLDS))
	 (force (getf coeffs 'NAVIER-STOKES::FORCE)))
    ;; loop over quadrature points
    (loop for shape-vals in (ip-values vecfe qrule) ; nr-comps x (n-basis x 1)
	  and shape-grads in (ip-gradients vecfe qrule) ; nr-comps x (n-basis x dim)
	  and global in (getf fe-geometry :global-coords)
	  and Dphi^-1 in (getf fe-geometry :gradient-inverses)
	  and weight in (getf fe-geometry :weights)
	  do
	  (let* ((sol-ip (and local-sol
			      (map 'vector #'(lambda (sv ls) (m* (transpose sv) ls))
				   local-sol shape-vals)))
		 (coeff-input (list :global global :solution sol-ip)))
	    (when viscosity (setq viscosity (evaluate viscosity coeff-input)))
	    (when reynolds (setq reynolds (evaluate reynolds coeff-input)))
	    (when force (setq force (evaluate force coeff-input)))
	    (when (= cell-dim dim)
	      (let* ((gradients (map 'vector (rcurry #'m* Dphi^-1) shape-grads)) ; nr-comps x (n-basis x dim)
		     (grad-p (vector-last gradients)))
		;; *** matrix ***
		(dotimes (i nr-comps)
		  (dotimes (j nr-comps)
		    ;; matrix
		    (when local-mat
		      ;; momentum part: tested with velocity
		      (when (< i dim)
			;; - \Delta u
			(when (= i j)
			  (unless (zerop viscosity)
			    (gemm! (* weight viscosity) (aref gradients i) (aref gradients j)
				   1.0d0 (aref local-mat i j) :NT)))
			;; (u \cdot \nabla) u
			(when (< j dim)
			  (unless (zerop reynolds)
			    (let* ((u-ip (apply #'matlisp::join-matrix
						(loop for u-comp across sol-ip
						      and i below dim collect u-comp)))
				   (u.grad_ui (m* (aref gradients i) u-ip)))
			      (gemm! (* weight reynolds)
				     u.grad_ui (aref shape-vals j)
				     1.0d0 (aref local-mat i j) :NT))))
			;; + \nabla p \vphi = - (\Div \phi) p
			(when (= j dim)
			  (gemm! weight (aref shape-vals i) (matrix-slice grad-p :from-col i :ncols 1)
				 1.0d0 (aref local-mat i j) :NT)
			  #+(or)
			  (gemm! (- weight) (matrix-slice (aref gradients i) :from-col i :ncols 1)
				 (aref shape-vals j) 1.0d0 (aref local-mat i j) :NT)
			  ))
		      ;; continuity part: tested with pressure
		      (when (= i dim)
			;; \Div u
			(when (< j dim)
			  (gemm! weight (aref shape-vals i)
				 (matrix-slice (aref gradients j) :from-col j :ncols 1)
				 1.0d0 (aref local-mat i j) :NT))))))))
	    ;; *** rhs ***
	    (dotimes (i nr-comps)
	      (when (and force local-rhs)
		(gemm! weight (aref shape-vals i) (aref force i)
		       1.0 (aref local-rhs i))))))))

	       

;;; Constraint assembly -> see system-fe.lisp

;;; Testing
(defun navier-stokes-fe-tests ()

  (let* ((level 1)
	 (problem (driven-cavity 2))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc  (navier-stokes-lagrange-fe 1 2 1)))
    (multiple-value-bind (matrix rhs)
	(discretize-globally problem h-mesh fedisc)
      (show rhs)
      (display matrix)))
  
  )

(tests::adjoin-femlisp-test 'navier-stokes-fe-tests)
