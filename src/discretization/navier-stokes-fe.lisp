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
(defpackage "FL.NAVIER-STOKES-FE"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM" "FL.NAVIER-STOKES" "FL.DISCRETIZATION")
  (:export "NAVIER-STOKES-LAGRANGE-FE")
  (:documentation "This package specializes the finite element
discretization for Navier-Stokes problems.  Up to now, we use only
generalized Taylor hood elements."))

(in-package "FL.NAVIER-STOKES-FE")

(defun navier-stokes-lagrange-fe (order dim delta)
  "Returns a generalized Taylor-Hood element @math{(Q^{k+delta})^d/Q^k} of
order @math{k} in dimension @math{d}."
  (lagrange-fe (map 'vector
		    (lambda (k)
		      (if (< k dim)
			  (+ order delta)
			  order))
		    (range<= 0 dim))))

#+(or)  ; is provided by ellsys discretization
(defvar *full-newton* t
  "If T, a full Newton linearization is used.  If NIL, the reaction term is
dropped.")

#+(or)  ; is provided by ellsys discretization
(defmethod discretize-locally ((problem <navier-stokes-problem>) coeffs vecfe qrule fe-geometry
			       &key matrix rhs local-u local-v
			       fe-parameters &allow-other-keys)
  "Local discretization for a Navier-Stokes problem."
  (assert (and (null local-u) (null local-v)))
  (let* ((nr-comps (nr-of-components vecfe))
	 (dim (1- nr-comps))
	 (cell (getf fe-geometry :cell))
	 (cell-dim (dimension cell))
	 (viscosity (get-coefficient coeffs 'FL.NAVIER-STOKES::VISCOSITY))
	 (reynolds (get-coefficient coeffs 'FL.NAVIER-STOKES::REYNOLDS))
	 (force (get-coefficient coeffs 'FL.NAVIER-STOKES::FORCE))
	 (solution (getf fe-parameters :solution)))
    ;; loop over quadrature points
    (loop
     for k from 0
     for shape-vals across (ip-values vecfe qrule) ; nr-comps x (n-basis x 1)
     and shape-grads across (ip-gradients vecfe qrule) ; nr-comps x (n-basis x dim)
     and global across (getf fe-geometry :global-coords)
     and Dphi across (getf fe-geometry :gradients)
     and Dphi^-1 across (getf fe-geometry :gradient-inverses)
     and weight across (getf fe-geometry :weights)
     do
     (let* ((coeff-input (construct-coeff-input
			  cell global Dphi shape-vals nil fe-parameters))
	    (sol-ip (getf coeff-input :solution)))
       (whereas ((force (and force (evaluate force coeff-input))))
	 ;; makes sense also for lower-dimensional cells
	 (when rhs
	   (dotimes (i nr-comps)
	     (gemm! weight (aref shape-vals i) (aref force i)
		    1.0 (aref rhs i)))))
       ;; the following should make sense only for dim-dimensional cells
       (when (= cell-dim dim)
	 (let* ((gradients (map 'vector (rcurry #'m* Dphi^-1) shape-grads)) ; nr-comps x (n-basis x dim)
		(grad-p (vector-last gradients)))
	   (whereas ((viscosity (and viscosity (evaluate viscosity coeff-input))))
	     ;; Stokes part: matrix 
	     (dotimes (i nr-comps)
	       (dotimes (j nr-comps)
		 ;; matrix
		 (when matrix
		   ;; momentum part: tested with velocity
		   (when (< i dim)
		     ;; - \Delta u
		     (when (= i j)
		       (when (and viscosity (not (zerop viscosity)))
			 (gemm! (* weight viscosity) (aref gradients i) (aref gradients j)
				1.0 (aref matrix i j) :NT)))
		     ;; + \phi \nabla p = - (\Div \phi) p
		     (when (= j dim)
		       (gemm! weight (aref shape-vals i) (matrix-slice grad-p :from-col i :ncols 1)
			      1.0 (aref matrix i j) :NT)))
		   ;; continuity part: tested with pressure
		   (when (= i dim)
		     ;; p \Div u
		     (when (< j dim)
		       (gemm! weight (aref shape-vals i)
			      (matrix-slice (aref gradients j) :from-col j :ncols 1)
			      1.0 (aref matrix i j) :NT)))))))
	   ;; Convection part
	   (whereas ((reynolds (and reynolds (evaluate reynolds coeff-input))))
	     (when (not (zerop reynolds))
	       (dbg :ns "Assembling Reynolds term")
	       (let ((u-ip (make-real-matrix
			    (loop for k below dim
				  collecting (vref (aref sol-ip k) 0)))))
		 (dotimes (i dim)
		   ;; (usol \cdot \nabla) u
		   (let ((u.grad_ui (m* (aref gradients i) u-ip)))
		     (gemm! (* weight reynolds)
			    (aref shape-vals i) u.grad_ui
			    1.0 (aref matrix i i) :NT))
		   ;; Full Newton linearization includes the term u (\nabla usol)
		   (when *full-newton*
		     (let ((grad-i (m*-tn (aref solution i) (aref gradients i))))
		       (dotimes (j dim)
			 (let ((factor (* weight reynolds (vref grad-i j))))
			   ;; put this term into the matrix
			   (gemm! factor (aref shape-vals i) (aref shape-vals j)
				  1.0 (aref matrix i j) :NT)
			   ;; and correct also the rhs
			   (when rhs
			     (axpy! (* factor (vref u-ip j)) (aref shape-vals i)
				    (aref rhs i)))))))))))
	   ))))))

;;; Constraint assembly -> see system-fe.lisp
	       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; GPS interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-discretization ((problem <navier-stokes-problem>) blackboard)
  (declare (ignore blackboard))
  (let* ((dim (dimension (domain problem)))
	 (order (or *suggested-discretization-order*
		    (if (<= dim 2) 2 1))))
    (navier-stokes-lagrange-fe order dim 1)))


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
    (let* ((dim 2) (domain (n-cube-domain dim))
	 (mesh (uniformly-refined-hierarchical-mesh domain 1))
	 (fe-class (navier-stokes-lagrange-fe 2 1 1))
	 (ansatz-space (make-instance '<ansatz-space> :fe-class fe-class :mesh mesh))
	 (I (interpolation-matrix ansatz-space))
	 (P (projection-matrix ansatz-space)))
    (assert (midentity-p (sparse-m* P I) 1.0e-10)))
  )

;;; (navier-stokes-fe-tests)
(fl.tests:adjoin-test 'navier-stokes-fe-tests)

