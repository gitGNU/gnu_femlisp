;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ellsys.lisp - Discretization of general elliptic systems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2006 Nicolas Neuss, University of Karlsruhe.
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
;;; NO EVENT SHALL THE AUTHOR, THE UNIVERSITY OF KARLSRUHE OR OTHER
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
(defpackage "FL.ELLSYS-FE"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM" "FL.ELLSYS"
	"FL.DISCRETIZATION")
  (:export)
  (:documentation "Finite element discretization of an elliptic system, see
the description in the ELLSYS package.  The result is a local matrix A and
local rhs b.  They will usually depend on u_old which is stored in the
solution vector.

This discretization is a unification of scalar and elasticity
discretization, and might replace those later on."))

(in-package :fl.ellsys-fe)

(defmethod discretize-locally
    ((problem <ellsys-problem>) coeffs vecfe qrule fe-geometry
     &key matrix rhs local-u local-v fe-parameters &allow-other-keys)
  "Local discretization for a pde system as described in the ELLSYS package
documentation."
  ;; and this situation should be checked before use
  (assert (and (null local-u) (null local-v)))

  (let ((diffusion-function (get-coefficient coeffs 'FL.ELLSYS::A))
	(convection-function (get-coefficient coeffs 'FL.ELLSYS::B))
	(reaction-function (get-coefficient coeffs 'FL.ELLSYS::C))
	(source-function (get-coefficient coeffs 'FL.ELLSYS::F))
	(g-source-function (get-coefficient coeffs 'FL.ELLSYS::G))
	(gamma-function (get-coefficient coeffs 'FL.ELLSYS::H))
	(cell (getf fe-geometry :cell))
	(ip-values (ip-values vecfe qrule))
	(ip-gradients (ip-gradients vecfe qrule)))

    ;; loop over quadrature points
    (loop
     for i from 0
     and global across (getf fe-geometry :global-coords)
     and Dphi across (getf fe-geometry :gradients)
     and Dphi^-1 across (getf fe-geometry :gradient-inverses)
     and weight across (getf fe-geometry :weights)
     do
     (let* ((shape-vals (aref ip-values i)) ; nr-comps x (n-basis x 1)
	    (shape-grads (aref ip-gradients i)) ; nr-comps x (n-basis x dim)
	    (gradients
	     (and Dphi^-1 (map 'vector (rcurry #'m* Dphi^-1) shape-grads)))
	    (coeff-input (construct-coeff-input
			  cell global Dphi shape-vals gradients fe-parameters))
	    (right-vals shape-vals)
	    (left-vals shape-vals)
	    (right-gradients gradients)
	    (left-gradients gradients))
       (let ((diff-tensor
		  (and diffusion-function
		       (evaluate diffusion-function coeff-input)))
	     (gamma (and gamma-function (evaluate gamma-function coeff-input))))
	 (when diff-tensor
	   ;; diffusion vector
	   (dbg :disc "Discretizing diffusion")
	   (for-each-key-and-entry
	    #'(lambda (indices diffusion)
		(destructuring-bind (i j) indices
		  (when (and diffusion (not (mzerop diffusion)))
		    (let ((fluxes (make-analog (aref gradients j))))
		      (gemm! 1.0 (aref left-gradients j) diffusion 0.0 fluxes)
		      (when matrix
			(gemm! weight fluxes (aref right-gradients i) 1.0 (mref matrix i j) :NT)
			;; gamma
			(when (and gamma rhs)
			  (gemm! weight fluxes (aref gamma j) 1.0 (aref rhs i))))))))
	    diff-tensor)))
       ;; convection
       (whereas ((velocity-tensor
		  (and convection-function
		       (evaluate convection-function coeff-input))))
	 (dbg :disc "Discretizing velocity")
	 (for-each-key-and-entry
	    #'(lambda (indices velocity)
		(destructuring-bind (i j) indices
		  (gemm! weight (m* (aref left-gradients i) velocity)
			 (aref right-vals j) 1.0 (mref matrix i j) :NT)))
	  velocity-tensor))
       ;; reaction
       (whereas ((reaction-tensor
		  (and reaction-function matrix
		       (evaluate reaction-function coeff-input))))
	 (dbg :disc "Discretizing reaction")
	 (for-each-key-and-entry
	  #'(lambda (indices reaction)
	      (destructuring-bind (i j) indices
		(unless (mzerop reaction)
		  (gemm! weight (m* (aref left-vals i) reaction) (aref right-vals j)
			 1.0 (mref matrix i j) :NT))))
	  reaction-tensor))
       ;; rhs-part / source
       (whereas ((source-vector
		  (and source-function rhs
		       (evaluate source-function coeff-input))))
	 (dbg :disc "Discretizing source")
	 (for-each-key-and-entry
	  #'(lambda (k source)
	      (gemm! weight (aref left-vals k) source
		     1.0 (aref rhs k)))
	  source-vector))
       (whereas ((g-source-vector
		  (and g-source-function rhs
		       (evaluate g-source-function coeff-input))))
	 (dbg :disc "Discretizing g-source")
	 (for-each-key-and-entry
	  #'(lambda (k g-source)
	      (gemm! weight (aref left-gradients k) g-source
		     1.0 (aref rhs k)))
	  g-source-vector))
       ))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; GPS interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-discretization ((problem <ellsys-problem>) blackboard)
  (declare (ignore blackboard))
  (let ((dim (dimension (domain problem))))
    (lagrange-fe (or *suggested-discretization-order*
		     (if (<= dim 2) 4 3))
		 :nr-comps (nr-of-components problem))))


;;; Testing

(defun discretize-ellsys (dim nr-comps order level)
  "Discretizes an elliptic system with @arg{nr-comps} components on the
@arg{dim}-dimensional unit cube with Lagrange fe of @arg{order}."
  (let* ((problem (ellsys-model-problem
		   dim nr-comps
		   :a (isotropic-diffusion dim #(1.0))
		   :f (coerce (loop repeat nr-comps collect #m(1.0)) 'vector)
		   :dirichlet (constraint-coefficient nr-comps 1)))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order :nr-comps nr-comps)))
    (discretize-globally problem h-mesh fedisc)))

(defun system-fe-tests ()
  
  (multiple-value-bind (matrix rhs)
      (discretize-ellsys 1 1 1 2)
    (fl.algebra:show matrix)
    (fl.algebra:show rhs))
  
  (time (discretize-ellsys 2 1 1 5))

  #+(or)
  (multiple-value-bind (matrix rhs)
      (discretize-ellsys 2 1 1 5)
    (fl.plot:plot (gesv matrix rhs) :component 0))

  #+(or)
  (multiple-value-bind (matrix rhs)
      (time (discretize-elasticity 2 1 5))
    (fl.plot:plot (gesv matrix rhs) :component 0))

  )

(fl.tests:adjoin-test 'system-fe-tests)




