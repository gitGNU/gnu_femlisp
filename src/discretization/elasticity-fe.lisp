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

(in-package :cl-user)
(defpackage "ELASTICITY-FE"
  (:use "COMMON-LISP" "FEMLISP.MATLISP" "MACROS" "UTILITIES" "ALGEBRA"
	"MESH" "PROBLEM" "ELASTICITY" "DISCRETIZATION")
  (:export))
(in-package :elasticity-fe)

(defmethod discretize-locally ((problem <elasticity-problem>) coeffs vecfe qrule fe-geometry
			       &key local-mat local-rhs local-sol local-u local-v)
  "Local discretization for an elasticity problem, in fact, even more
general it can be used for solving elliptic systems.  The basic idea is
that we work with arrays having standard-matrix entries whereas in the
scalar case we worked directly with standard-matrix."
  (declare (type (array real-matrix (* *)) local-mat)
	   (type (or null (array real-matrix (*)))
		 local-sol local-rhs local-u local-v))
  ;; we want isotropic ansatz spaces...
  (assert (same-p (components vecfe)))
  ;; and this situation should be checked before use
  (assert (and (null local-u) (null local-v)))
  
  (let ((fe (aref (components vecfe) 0))
	(nr-comps (nr-of-components vecfe))
	(elasticity-tensor (getf coeffs 'ELASTICITY::ELASTICITY))
	(gamma-function (getf coeffs 'ELASTICITY::GAMMA))
	(reaction-function (getf coeffs 'ELASTICITY::REACTION))
	(force-function (getf coeffs 'ELASTICITY::FORCE)))

    ;; loop over quadrature points
    (loop for shape-vals in (ip-values fe qrule) ; (n-basis x 1)-matrix
	  and shape-grads in (ip-gradients fe qrule) ; (n-basis x dim)-matrix
	  and global in (getf fe-geometry :global-coords)
	  and Dphi^-1 in (getf fe-geometry :gradient-inverses)
	  and weight in (getf fe-geometry :weights)
	  do
	  (let* ((transposed-sv (transpose shape-vals))
		 (sol-ip (and local-sol (map 'vector #'(lambda (ls) (m* transposed-sv ls))
					     local-sol)))
		 (coeff-input (list :global global :solution sol-ip))
		 (ip-tensor (and elasticity-tensor
				 (evaluate elasticity-tensor coeff-input)))
		 (gamma (and gamma-function local-rhs
			     (evaluate gamma-function coeff-input)))
		 (reaction (and reaction-function
				(evaluate reaction-function coeff-input)))
		 (force (and force-function (evaluate force-function coeff-input)))
		 (gradients (m* shape-grads Dphi^-1)) ; (n-basis x dim)-matrix
		 (right-gradients gradients)
		 (left-gradients gradients)
		 (fluxes (make-analog gradients))
		 (right-vals shape-vals)
		 (left-vals shape-vals))
	    
	    (dotimes (i nr-comps)
	      ;; matrix-part
	      (dotimes (j nr-comps)
		;; diffusion 
		(when ip-tensor
		  (let ((D (aref ip-tensor i j)))
		    (unless (mzerop D)
		      (gemm! 1.0 left-gradients D 0.0 fluxes)
		      (when local-mat
			(gemm! weight fluxes right-gradients 1.0d0 (aref local-mat i j) :NT))
		      ;; gamma
		      (when (and gamma local-rhs)
			(gemm! weight fluxes (aref gamma j) 1.0d0 (aref local-rhs i))))))
		;; reaction
		(when (and reaction local-mat)
		  (let ((R (aref reaction i j)))
		    (unless (mzerop R)
		      (gemm! (* weight R)
			     left-vals right-vals 1.0d0 local-mat :NT)))))
	      ;; rhs-part / force
	      (when (and force local-rhs)
		(m+! (m* left-vals (scal weight (aref force i)))
		     (aref local-rhs i))))))))

;;; Assembly of boundary conditions -> system-fe.lisp

;;; Testing
(defun elasticity-fe-tests ()
  (let* ((dim 1) (order 1) (level 2)
	 (problem (standard-elasticity-problem
		   (n-cube-domain dim) :lambda 1.0 :mu 1.0
		   :force (constantly (unit-vector dim 0))))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order :nr-comps dim)))
    (multiple-value-bind (mat rhs)
	(discretize-globally problem h-mesh fedisc)
      (m* (sparse-ldu mat) rhs)))
  )

(tests::adjoin-femlisp-test 'elasticity-fe-tests)
