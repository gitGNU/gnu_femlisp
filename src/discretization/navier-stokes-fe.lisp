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
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES"
	"FL.MATLISP" "ALGEBRA" "FL.FUNCTION"
	"MESH" "PROBLEM" "NAVIER-STOKES" "DISCRETIZATION" "NAVIER-STOKES")
  (:export "NAVIER-STOKES-LAGRANGE-FE"))
(in-package :navier-stokes-fe)

(defclass <navier-stokes-lagrange-fe> (<vector-fe-discretization>)
  ()
  (:documentation "Up to now these are the generalized Taylor-Hood elements
(Q_{k+1})^dim/Q_k."))

(defun navier-stokes-lagrange-fe (order dim delta)
  "Returns (Q_{k+delta})^dim/Q_k finite elements."
  (make-instance
   '<navier-stokes-lagrange-fe>
   :components
   (append
    (make-list dim :initial-element (lagrange-fe (+ order delta)))
    (list (lagrange-fe order :integration-order (+ order delta))))))

(defvar *full-newton* t
  "If T, a full Newton linearization is used.  If NIL, the reaction term is
dropped.")

(defmethod discretize-locally ((problem <navier-stokes-problem>) coeffs vecfe qrule fe-geometry
			       &key local-mat local-rhs local-sol local-u local-v
			       coefficient-parameters &allow-other-keys)
  "Local discretization for a Navier-Stokes problem."
  (assert (and (null local-u) (null local-v)))
  
  (let* ((nr-comps (nr-of-components vecfe))
	 (dim (1- nr-comps))
	 (cell-dim (dimension (getf fe-geometry :cell)))
	 (viscosity (getf coeffs 'NAVIER-STOKES::VISCOSITY))
	 (reynolds (getf coeffs 'NAVIER-STOKES::REYNOLDS))
	 (force (getf coeffs 'NAVIER-STOKES::FORCE)))
    ;; loop over quadrature points
    (loop
     for k from 0
     for shape-vals across (ip-values vecfe qrule) ; nr-comps x (n-basis x 1)
     and shape-grads across (ip-gradients vecfe qrule) ; nr-comps x (n-basis x dim)
     and global in (getf fe-geometry :global-coords)
     and Dphi^-1 in (getf fe-geometry :gradient-inverses)
     and weight in (getf fe-geometry :weights)
	  do
	  (let* ((sol-ip (and local-sol (map 'vector #'m*-tn local-sol shape-vals)))
		 (coeff-input
		  (list* :global global :solution sol-ip
			 (loop for (key data) on coefficient-parameters by #'cddr
			       collect key collect (aref data k)))))
	    (when viscosity (setq viscosity (evaluate viscosity coeff-input)))
	    (when reynolds (setq reynolds (evaluate reynolds coeff-input)))
	    (when force (setq force (evaluate force coeff-input)))
	    (when (= cell-dim dim)
	      (let* ((gradients (map 'vector (rcurry #'m* Dphi^-1) shape-grads)) ; nr-comps x (n-basis x dim)
		     (grad-p (vector-last gradients)))
		;; Stokes part: matrix 
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
				   1.0 (aref local-mat i j) :NT)))
			;; + \phi \nabla p = - (\Div \phi) p
			(when (= j dim)
			  (gemm! weight (aref shape-vals i) (matrix-slice grad-p :from-col i :ncols 1)
				 1.0 (aref local-mat i j) :NT)))
		      ;; continuity part: tested with pressure
		      (when (= i dim)
			;; p \Div u
			(when (< j dim)
			  (gemm! weight (aref shape-vals i)
				 (matrix-slice (aref gradients j) :from-col j :ncols 1)
				 1.0 (aref local-mat i j) :NT))))))
		;; Stokes part: rhs
		(dotimes (i nr-comps)
		  (when (and force local-rhs)
		    (gemm! weight (aref shape-vals i) (aref force i)
			   1.0 (aref local-rhs i))))
		;; Convection part
		(unless (zerop reynolds)
		  (let ((u-ip (make-real-matrix
			       (loop for k below dim
				     collecting (vref (aref sol-ip k) 0)))))
		    (dotimes (i dim)
		      ;; (usol \cdot \nabla) u
		      (let ((u.grad_ui (m* (aref gradients i) u-ip)))
			(gemm! (* weight reynolds)
			       (aref shape-vals i) u.grad_ui
			       1.0 (aref local-mat i i) :NT))
		      ;; Full Newton linearization includes the term u (\nabla usol)
		      (when *full-newton*
			(let ((grad-i (m*-tn (aref local-sol i) (aref gradients i))))
			  (dotimes (j dim)
			    (let ((factor (* weight reynolds (vref grad-i j))))
			      ;; put this term into the matrix
			      (gemm! factor (aref shape-vals i) (aref shape-vals j)
				     1.0 (aref local-mat i j) :NT)
			      ;; and correct also the rhs
			      (when local-rhs
				(axpy! (* factor (vref u-ip j)) (aref shape-vals i)
				       (aref local-rhs i))))))))))
		))))))

;;; Constraint assembly -> see system-fe.lisp
	       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; GPS interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-discretization ((problem <navier-stokes-problem>) blackboard)
  (let* ((dim (dimension (domain problem)))
	 (order (if (<= dim 2) 2 1)))
    (navier-stokes-fe::navier-stokes-lagrange-fe order dim 1)))


;;; Testing
(defun navier-stokes-fe-tests ()
  (let* ((level 1)
	 (problem (driven-cavity 2))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (navier-stokes-lagrange-fe 1 2 1)))
    (multiple-value-bind (matrix rhs)
	(discretize-globally problem h-mesh fedisc)
      (show rhs)
      (display matrix)))
  )

(fl.tests:adjoin-test 'navier-stokes-fe-tests)
