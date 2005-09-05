;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  cdr-fe.lisp
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
(defpackage "FL.CDR-FE"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.DEBUG"
	"FL.MATLISP" "FL.ALGEBRA" "FL.FUNCTION"
	"FL.MESH" "FL.PROBLEM" "FL.CDR" "FL.DISCRETIZATION")
  (:export )
  (:documentation "This package specializes the finite element
discretization for convection-diffusion-reaction problems.

It can handle the following equation:

@math{- \partial_i (K_{ij} (\partial_j u + g_j) + \partial_i (c_i u) + r u = f}

Here, @math{K} is the diffusion tensor, @math{c} is the convection vector,
@math{r} is the reaction coefficient, @math{f} is the source function, and
@math{g} is a distributional source."))

(in-package "FL.CDR-FE")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Local assembly
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod discretize-locally ((problem <cdr-problem>) coeffs fe qrule fe-geometry
			       &key local-mat local-rhs local-sol local-u local-v
			       mass-factor stiffness-factor
			       coefficient-parameters)
  "Local discretization for a convection-diffusion-reaction equation."
  ;; extract active coefficient functions
  (let ((diffusion (getf coeffs 'FL.CDR::DIFFUSION))
	(convection (getf coeffs 'FL.CDR::CONVECTION))
	(gamma (getf coeffs 'FL.CDR::GAMMA))
	(source (getf coeffs 'FL.CDR::SOURCE))
	(reaction (getf coeffs 'FL.CDR::REACTION))
	(cell (getf fe-geometry :cell)))
    (dbg :fe "Local discretization on cell ~A" cell)
    
    ;; loop over quadrature points
    (loop
     for i from 0
     and shape-vals across (ip-values fe qrule) ; (n-basis x 1)-matrix
     and shape-grads across (ip-gradients fe qrule) ; (n-basis x dim)-matrix
     and global in (getf fe-geometry :global-coords)
     and Dphi in (getf fe-geometry :gradients)
     and Dphi^-1 in (getf fe-geometry :gradient-inverses)
     and weight in (getf fe-geometry :weights)
     do
     (let* ((right-vals (if local-u (dot local-u shape-vals) shape-vals))
	    (left-vals (if local-v (dot local-v shape-vals) shape-vals))
	    (sol-ip (and local-sol (m*-tn shape-vals local-sol)))
	    (coeff-input
	     (list* :global global :solution sol-ip :Dphi Dphi :cell cell
		    (loop for (key data) on coefficient-parameters by #'cddr
			  collect key collect (aref data i)))))
       
       (when (or (and local-mat (or diffusion convection)
		      (not (zerop stiffness-factor)))
		 (and local-rhs diffusion gamma))
	 (let* ((gradients (m* shape-grads Dphi^-1)) ; (n-basis x dim)-matrix
		(right-gradients (if local-u (m*-tn local-u gradients) gradients))
		(left-gradients (if local-v (m*-tn local-v gradients) gradients)))
	   ;; diffusion
	   (when diffusion
	     (let* ((diff-tensor (ensure-matlisp (evaluate diffusion coeff-input)))
		    (fluxes (m* left-gradients diff-tensor))) ; (n-basis x dim)-matrix
	       (when local-mat
		 (gemm! weight fluxes right-gradients 1.0 local-mat :NT))
	       ;; gamma
	       (when (and gamma local-rhs)
		 (gemm! weight fluxes (evaluate gamma coeff-input)
			1.0 local-rhs))
	       ))
	   ;; convection
	   (when (and convection local-mat)
	     (let* ((velocity-vector (evaluate convection coeff-input)))
	       (gemm! (- weight) (m* left-gradients velocity-vector) right-vals 1.0 local-mat :NT)))
	   ))
	    
       ;; reaction
       (when (and local-mat reaction (not (zerop stiffness-factor)))
	 (let* ((reaction (evaluate reaction coeff-input))
		(factor (* weight reaction)))
	   (gemm! factor left-vals right-vals 1.0 local-mat :NT)))
       
       ;; mass matrix
       (when (and local-mat (not (zerop mass-factor)))
	 (let* ((factor (* weight mass-factor)))
	   (gemm! factor left-vals right-vals 1.0 local-mat :NT)))
       
       ;; source
       (when (and source local-rhs)
	 (let ((source (evaluate source coeff-input)))
	   (when (numberp source) (setq source (make-real-matrix `((,source)))))
	   (gemm! weight left-vals source 1.0 local-rhs)))))
    
    ;; custom fe rhs
    (whereas ((fe-rhs (and local-rhs (getf coeffs 'FL.CDR::FE-RHS))))
      (m+! (evaluate fe-rhs (list :cell cell :fe fe))
	   local-rhs))
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Assemble (Dirichlet) boundary
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod essential-boundary-constraints ((problem <cdr-problem>) (ansatz-space <ansatz-space>)
					    &key level (where :surface) interface)
  (assert (or (not (eq where :surface)) interface))
  (let* ((problem (problem ansatz-space))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (level-skel (if level (cells-on-level h-mesh level) h-mesh))
	 (fe-class (fe-class ansatz-space))
	 (constraints-P (make-ansatz-space-automorphism ansatz-space))
	 (constraints-Q (make-ansatz-space-automorphism ansatz-space))
	 (constraints-rhs (make-ansatz-space-vector ansatz-space))
	 (multiplicity (multiplicity ansatz-space)))
    (doskel (cell level-skel)
      (when (ecase where
	      (:refined (refined-p cell h-mesh))
	      (:surface (or (not (refined-p cell h-mesh))
			    (member-of-skeleton? cell interface)))
	      (:all t))
	(let* ((coeffs (coefficients-of-cell cell h-mesh problem))
	       (dirichlet-function (getf coeffs 'FL.CDR::CONSTRAINT))
	       (cell-key (cell-key cell h-mesh)))
	  (when dirichlet-function
	    (let ((fe (get-fe fe-class cell)))
	      (loop+ ((dof (fe-dofs fe))
		      (j (range :below (nr-of-inner-dofs fe))))
		 do
		 (let ((k (dof-in-vblock-index dof)))
		   ;; The following is only correct for degrees of freedom of
		   ;; Lagrange type.  Perhaps one should use Hermite finite
		   ;; cells only in the interior?
		   (setf (mref (mref constraints-P cell-key cell-key) k k) 1.0)
		   ;;		       (clear-row mat cell k)
		   ;;		       (setf (mref (mref mat cell cell) k k) 1.0)
		   (let* ((value
			   (evaluate dirichlet-function
				     (list :local (dof-coord dof)
					   :global (local->global cell (dof-gcoord dof)))))
			  (values (cond ((numberp value)
					 (scal value (ones 1 multiplicity)))
					((vectorp value) (aref value 0))
					(t value))))
		     (minject values (vref constraints-rhs cell-key)
			      k 0)))))))))
    (values constraints-P constraints-Q constraints-rhs)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; GPS interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-discretization ((problem <cdr-problem>) blackboard)
  (declare (ignore blackboard))
  (let* ((dim (dimension (domain problem)))
	 (order (or *suggested-discretization-order*
		    (cond ((<= dim 2) 4)
			  ((= dim 3) 3)
			  (t 1)))))
    (lagrange-fe order)))


;;; Testing
(defun cdr-fe-tests ()
  (let* ((order 1) (level 1)
	 (problem (cdr-model-problem 1))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order))
	 (as (make-fe-ansatz-space fedisc problem h-mesh)))
    (with-items (&key matrix rhs)
	(fe-discretize
	 (blackboard :ansatz-space as
		     :mass-factor 1.0 :stiffness-factor 1.0))
      (show (getrs (sparse-ldu matrix) rhs))))
  
  (let* ((level 1) (order 1)
	 (problem (cdr-model-problem 1))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (as1 (make-fe-ansatz-space (lagrange-fe order) problem h-mesh))
	 (as2 (make-fe-ansatz-space (lagrange-fe (1+ order)) problem h-mesh)))
    (with-items (&key matrix rhs)
	(fe-discretize (blackboard :ansatz-space as1))
      ;; set constraints
      (fe-discretize (blackboard :ansatz-space as2))
      (let* ((sol (getrs (sparse-ldu matrix) rhs))
	     (low->high (transfer-matrix as1 as2))
	     (high->low (transfer-matrix as2 as1)))
	(assert (< (norm (m- sol (m* high->low (m* low->high sol)))) 1.0e-10)))))
  )

(fl.tests:adjoin-test 'cdr-fe-tests)

