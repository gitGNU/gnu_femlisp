;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  cdrsys-fe.lisp
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

(in-package "COMMON-LISP-USER")
(defpackage "FL.CDRSYS-FE"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM" "FL.CDRSYS" "FL.DISCRETIZATION")
  (:export)
  (:documentation "This package specializes the finite element
discretization for systems of convection-diffusion-reaction type."))

(in-package "FL.CDRSYS-FE")

(defmethod discretize-locally
    (<cdrsys-problem> coeffs vecfe qrule fe-geometry
     &key matrix rhs local-u local-v fe-parameters &allow-other-keys)
  "Local discretization for a convection-diffusion-reaction system.
Derived from the elasticity discretization.  A unification is intended for
the future."
  (declare (optimize debug))
  ;; we need isotropic ansatz spaces...
  (assert (same-p (components vecfe)))
  ;; and this situation should be checked before use
  (assert (and (null local-u) (null local-v)))
  
  (let ((diffusion-function (get-coefficient coeffs 'FL.CDRSYS::DIFFUSION))
	(convection-function (get-coefficient coeffs 'FL.CDRSYS::CONVECTION))
	(source-function (get-coefficient coeffs 'FL.CDRSYS::SOURCE))
	(reaction-function (get-coefficient coeffs 'FL.CDRSYS::REACTION))
	(cell (getf fe-geometry :cell)))

    ;; loop over quadrature points
    (loop
     for shape-vals across (ip-values vecfe qrule) ; nr-comps x (n-basis x 1)
     and shape-grads across (ip-gradients vecfe qrule) ; nr-comps x (n-basis x dim)
     and global in (getf fe-geometry :global-coords)
     and Dphi in (getf fe-geometry :gradients)
     and Dphi^-1 in (getf fe-geometry :gradient-inverses)
     and weight in (getf fe-geometry :weights)
     do
     (let* ((gradients
	     (and Dphi^-1 (map 'vector (rcurry #'m* Dphi^-1) shape-grads)))
	    (coeff-input (construct-coeff-input
			  cell global Dphi shape-vals gradients fe-parameters))
	    (right-vals shape-vals)
	    (left-vals shape-vals)
	    (right-gradients gradients)
	    (left-gradients gradients))
       (whereas ((diff-tensor
		  (and diffusion-function
		       (evaluate diffusion-function coeff-input))))
	 ;; diffusion vector
	 (for-each-key-and-entry
	  #'(lambda (k diffusion)
	      (when (and diffusion (not (mzerop diffusion)))
		(let ((fluxes (make-analog (aref gradients k))))
		  (gemm! 1.0 (aref left-gradients k) diffusion 0.0 fluxes)
		  (when matrix
		    (gemm! weight fluxes (aref right-gradients k) 1.0 (mref matrix k k) :NT)))))
	  diff-tensor))
       ;; convection
       (whereas ((velocity-vector
		  (and convection-function
		       (evaluate convection-function coeff-input))))
	 (for-each-key-and-entry
	  #'(lambda (k velocity)
	      (gemm! (- weight) (m* (aref left-gradients k) velocity)
		     (aref right-vals k) 1.0 (mref matrix k k) :NT))
	  velocity-vector))
       ;; reaction
       (whereas ((reaction-matrix
		  (and reaction-function matrix
		       (evaluate reaction-function coeff-input))))
	 (for-each-key-and-entry
	  #'(lambda (k reaction)
	      (unless (mzerop reaction)
		(gemm! (* weight reaction)
		       (aref left-vals k) (aref right-vals k)
		       1.0 (mref matrix k k) :NT)
		))
	  reaction-matrix))
       ;; rhs-part / source
       (whereas ((source-vector
		  (and source-function rhs
		       (evaluate source-function coeff-input))))
	 (for-each-key-and-entry
	  #'(lambda (k source)
	      (gemm! weight (aref left-vals k) source
		     1.0 (aref rhs k)))
	  source-vector))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; GPS interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-discretization ((problem <cdrsys-problem>) blackboard)
  (let ((dim (dimension (domain problem))))
    (lagrange-fe (or *suggested-discretization-order*
		     (if (<= dim 2) 4 3))
		 :nr-comps (nr-of-components problem))))

;;; Testing
(defun cdrsys-fe-tests ()
  #+(or) ; something like the following
  (let* ((order 1) (level 2)
	 (problem (cdr-model-problem 1))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (multiple-value-bind (matrix rhs)
	(discretize-globally problem h-mesh fedisc)
      (getrs (sparse-ldu matrix) rhs)))
  )

(fl.tests:adjoin-test 'cdrsys-fe-tests)




