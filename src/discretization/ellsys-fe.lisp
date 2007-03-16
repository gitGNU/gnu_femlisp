;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ellsys-fe.lisp - Discretization of general elliptic systems
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
  (:documentation "Finite element discretization of a general second order
elliptic system, see the description in the ELLSYS package.  The result is
a local matrix A and local rhs b.  They will usually depend on u_old which
is stored in the solution vector."))

(in-package :fl.ellsys-fe)

(defmethod discretize-locally
    ((problem <ellsys-problem>) coeffs vecfe qrule fe-geometry
     &key matrix rhs local-u local-v fe-parameters &allow-other-keys)
  "Local discretization for a pde system of the form described in the
documentation of the package @package{ELLSYS}."

  (assert (and (null local-u) (null local-v)))  ; NYI

  (let ((cell (getf fe-geometry :cell))
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
	    (gradients (and Dphi^-1 (map 'vector (rcurry #'m* Dphi^-1) shape-grads)))
	    (coeff-input (construct-coeff-input
			  cell global Dphi shape-vals gradients fe-parameters))
	    (right-vals shape-vals)
	    (left-vals shape-vals)
	    (right-gradients gradients)
	    (left-gradients gradients))
       (assert (vectorp left-vals))
       (loop
	for coeff in coeffs
	for name = (coefficient-name coeff) do
	;; matrix part
	(when matrix
	  (case name
	    (FL.ELLSYS::A  ; diffusion
	     (let ((diff-tensor (evaluate coeff coeff-input)))
	       (dbg :disc "Discretizing diffusion")
	       (for-each-entry-and-key
		#'(lambda (diffusion i j)
		    (gemm! weight (m* (aref left-gradients j) diffusion) (aref right-gradients i)
			   1.0 (mref matrix i j) :NT))
		diff-tensor)))
	    (FL.ELLSYS::B
	     (dbg :disc "Discretizing convection")
	     (let ((velocity-tensor (evaluate coeff coeff-input)))
	       (for-each-entry-and-key
		#'(lambda (velocity i j)
		    (gemm! weight (m* (aref left-gradients i) velocity)
			   (aref right-vals j) 1.0 (mref matrix i j) :NT))
		velocity-tensor)))
	    (FL.ELLSYS::C
	     (dbg :disc "Discretizing non-conservative convection")
	     (let ((velocity-tensor (evaluate coeff coeff-input)))
	       (for-each-entry-and-key
		#'(lambda (velocity i j)
		    (gemm! weight (aref left-vals i)
			   (m* (aref right-gradients j) velocity)
			   1.0 (mref matrix i j) :NT))
		velocity-tensor)))
	    (FL.ELLSYS::R
	     (dbg :disc "Discretizing reaction")
	     (let ((reaction-tensor (evaluate coeff coeff-input)))
	       (for-each-entry-and-key
		#'(lambda (reaction i j)
		    (unless (mzerop reaction)
		      (when (standard-matrix-p reaction)
			(assert (= 1 (nrows reaction) (ncols reaction)))
			(setq reaction (vref reaction 0)))
		      (gemm! (* weight reaction) (aref left-vals i) (aref right-vals j)
			     1.0 (mref matrix i j) :NT)))
		reaction-tensor)))))
	(when rhs
	  (case name
	    (FL.ELLSYS::F
	     (dbg :disc "Discretizing standard source")
	     (let ((source-vector  (evaluate coeff coeff-input)))
	       (for-each-entry-and-key
		#'(lambda (source k)
		    (gemm! weight (aref left-vals k) (ensure-matlisp source)
			   1.0 (aref rhs k)))
		source-vector)))
	    (FL.ELLSYS::G
	     (dbg :disc "Discretizing g-source")
	     (let ((g-source-vector (evaluate coeff coeff-input)))
	       (for-each-entry-and-key
		#'(lambda (g-source k)
		    (gemm! weight (aref left-gradients k) g-source
			   1.0 (aref rhs k)))
		g-source-vector)))
	    (FL.ELLSYS::H
	     (dbg :disc "Discretizing h-source")
	     (let ((diff-tensor (evaluate (get-coefficient coeffs 'FL.ELLSYS::A) coeff-input))
		   (gamma (evaluate coeff coeff-input)))
	       (for-each-entry-and-key
		#'(lambda (diffusion i j)
		    (gemm! weight (m* (aref left-gradients j) diffusion) (aref gamma j)
			   1.0 (aref rhs i)))
		diff-tensor))))))))
    ;; nonstandard rhs (outside of ip loop)
    (whereas ((fe-rhs-function (and rhs (get-coefficient coeffs 'FL.ELLSYS::FE-RHS))))
      (m+! (evaluate fe-rhs-function (list :cell cell :fe vecfe :geometry fe-geometry))
	   rhs))
    ))

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
		   dim `((u ,nr-comps))
		   :a (isotropic-diffusion dim #(1.0))
		   :f (coerce (loop repeat nr-comps collect #m(1.0)) 'vector)
		   :dirichlet (constraint-coefficient nr-comps 1)))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order :nr-comps nr-comps)))
    (discretize-globally problem h-mesh fedisc)))

;;;; Testing

(defun ellsys-fe-tests ()
  (dbg-on :disc)
  (multiple-value-bind (matrix rhs)
      (discretize-ellsys 1 1 1 2)
    (fl.algebra:show matrix)
    (fl.algebra:show rhs)
    (assert (zerop (norm (m- (fe-value (gesv matrix rhs) #(0.5))
			     (vector #m(0.125)))))))
  (dbg-off)

  (let* ((dim 1) (order 1) (level 2)
	 (problem (fl.cdr::cdr-model-problem dim :diffusion nil))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (discretize-globally problem h-mesh fedisc))
  
  (time (multiple-value-bind (matrix rhs)
	    (discretize-ellsys 2 1 1 5)
	  (fe-value (gesv matrix rhs) #d(0.5 0.5))))

  #+(or)  ; plot is not yet available when loading
  (multiple-value-bind (matrix rhs)
      (discretize-ellsys 2 1 1 5)
    (fl.plot:plot (gesv matrix rhs)))
  
  #+(or)
  (multiple-value-bind (matrix rhs)
      (time (discretize-elasticity 2 1 5))
    (fl.plot:plot (gesv matrix rhs) :component 0))

  )

;;; (ellsys-fe-tests)
(fl.tests:adjoin-test 'ellsys-fe-tests)
