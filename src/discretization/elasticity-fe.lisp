;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  elasticity.lisp
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

(in-package :CL-USER)
(defpackage "FL.ELASTICITY-FE"
  (:use "COMMON-LISP" "FL.MATLISP" "FL.MACROS" "FL.UTILITIES" "FL.ALGEBRA" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM" "FL.ELASTICITY" "FL.DISCRETIZATION")
  (:export)
  (:documentation "This package specializes the finite element
discretization for elasticity problems.  In fact, even more general
elliptic systems can be handled by providing suitable coefficient matrices.
Here, local solution and local right-hand side are vectors of standard
matrices (which allows multiple solutions), and the local stiffness matrix
is an array of rank 2 with standard-matrix entries."))
(in-package "FL.ELASTICITY-FE")

(defmethod discretize-locally
    ((problem <elasticity-problem>) coeffs vecfe qrule fe-geometry
     &key matrix rhs local-u local-v fe-parameters &allow-other-keys)
  "Local discretization for an elasticity problem."
  ;; we want isotropic ansatz spaces...
  (assert (same-p (components vecfe)))
  ;; and this situation should be checked before use
  (assert (and (null local-u) (null local-v)))
  
  (let ((fe (aref (components vecfe) 0))
	(nr-comps (nr-of-components vecfe))
	(elasticity-tensor (get-coefficient coeffs 'FL.ELASTICITY::ELASTICITY))
	(gamma-function (get-coefficient coeffs 'FL.ELASTICITY::GAMMA))
	(reaction-function (get-coefficient coeffs 'FL.ELASTICITY::REACTION))
	(force-function (get-coefficient coeffs 'FL.ELASTICITY::FORCE))
	(cell (getf fe-geometry :cell)))

    ;; loop over quadrature points
    (loop
     for k from 0
     and shape-vals across (ip-values fe qrule) ; (n-basis x 1)-matrix
     and shape-grads across (ip-gradients fe qrule) ; (n-basis x dim)-matrix
     and global in (getf fe-geometry :global-coords)
     and Dphi^-1 in (getf fe-geometry :gradient-inverses)
     and weight in (getf fe-geometry :weights)
     do
     (let* ((gradients (and Dphi^-1 (m* shape-grads Dphi^-1)))
	    (coeff-input (construct-coeff-input
			  cell global shape-vals gradients fe-parameters))
	    (ip-tensor (and elasticity-tensor
			    (evaluate elasticity-tensor coeff-input)))
	    (gamma (and gamma-function rhs
			(evaluate gamma-function coeff-input)))
	    (reaction (and reaction-function
			   (evaluate reaction-function coeff-input)))
	    (force (and force-function (evaluate force-function coeff-input)))
	    (right-vals shape-vals)
	    (left-vals shape-vals)
	    (right-gradients gradients)
	    (left-gradients gradients))
       (when ip-tensor
	 (dotimes (i nr-comps)
	   (dotimes (j nr-comps)
	     ;; diffusion 
	     (let ((D (aref ip-tensor i j)))
	       (unless (mzerop D)
		 (let ((fluxes (make-analog gradients)))
		   (gemm! 1.0 left-gradients D 0.0 fluxes)
		   (when matrix
		     (gemm! weight fluxes right-gradients 1.0 (aref matrix i j) :NT))
		   ;; gamma
		   (when (and gamma rhs)
		     (gemm! weight fluxes (aref gamma j) 1.0 (aref rhs i)))))))))
       ;; reaction
       (when (and reaction matrix)
	 (dotimes (i nr-comps)
	   (dotimes (j nr-comps)
	     (let ((R (aref reaction i j)))
	       (unless (mzerop R)
		 (gemm! (* weight R)
			left-vals right-vals 1.0 matrix :NT))))))
       ;; rhs-part / force
       (when (and force rhs)
	 (dotimes (i nr-comps)
	   (gemm! weight left-vals (aref force i)
		  1.0 (aref rhs i))))))))

;;; Assembly of boundary conditions -> system-fe.lisp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; GPS interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-discretization ((problem <elasticity-problem>) blackboard)
  (declare (ignore blackboard))
  (let* ((dim (dimension (domain problem)))
	 (order (or *suggested-discretization-order*
		    (if (<= dim 2) 4 3))))
    (lagrange-fe order :nr-comps dim)))


;;; Testing
(defun elasticity-fe-tests ()
  (let* ((dim 1) (order 1) (level 1)
	 (problem (standard-elasticity-problem
		   (n-cube-domain dim) :lambda 1.0 :mu 1.0
		   :force (constantly (vector #m(1.0)))))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order :nr-comps dim)))
    (multiple-value-bind (mat rhs)
	(discretize-globally problem h-mesh fedisc)
      (getrs (sparse-ldu mat) rhs) :component 0))
  )

;;; (elasticity-fe-tests)
(fl.tests:adjoin-test 'elasticity-fe-tests)
