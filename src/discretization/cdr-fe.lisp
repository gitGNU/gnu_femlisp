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
(defpackage "CDR-FE"
  (:use "COMMON-LISP" "FL.MACROS" "FL.UTILITIES" "FL.MATLISP"
	"MESH" "PROBLEM" "CDR" "DISCRETIZATION" "ALGEBRA" "FL.FUNCTION")
  (:export ))
(in-package :cdr-fe)

#|
Local and boundary assembly for finite element discretizations of
convection-diffusion-reaction problems.  It can handle the following
equation (D_i=\partial{x_i}):

- D_i (K_ij (D_j u + g_j) + D_i (c u) + alpha r u = alpha f

Here, K is the diffusion tensor, c is convection, r is reaction and f is
source.  alpha is a coefficient which arises from scaling of u in porous
media applications.  It is put separately from r and f, because it is used
also inside time-stepping schemes as multiplier of 1/delta_t.

Ideally, this module should be the only one which has to be changed for
other problems or other finite element discretizations.  |#

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Local assembly
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod discretize-locally ((problem <cdr-problem>) coeffs fe qrule fe-geometry
			       &key local-mat local-rhs local-sol local-u local-v)
  "Local discretization for a convection-diffusion-reaction equation."
  #+(or)
  (declare (type (or null real-matrix) local-mat)
	   (type (or null real-matrix) local-sol local-rhs local-u local-v))
  (let ((diffusion-function (getf coeffs 'CDR::DIFFUSION))
	(convection-function (getf coeffs 'CDR::CONVECTION))
	(gamma-function (getf coeffs 'CDR::GAMMA))
	(alpha-function (getf coeffs 'CDR::ALPHA))
	(source-function (getf coeffs 'CDR::SOURCE))
	(reaction-function (getf coeffs 'CDR::REACTION)))
    ;; loop over quadrature points
    (loop for shape-vals in (ip-values fe qrule) ; (n-basis x 1)-matrix
	  and shape-grads in (ip-gradients fe qrule) ; (n-basis x dim)-matrix
	  and global in (getf fe-geometry :global-coords)
	  and Dphi^-1 in (getf fe-geometry :gradient-inverses)
	  and weight in (getf fe-geometry :weights)
	  do
	  (let* ((transposed-sv (transpose shape-vals))
		 (sol-ip (and local-sol (m* transposed-sv local-sol)))
		 (right-vals (if local-u (dot local-u shape-vals) shape-vals))
		 (left-vals (if local-v (dot local-v shape-vals) shape-vals))
		 (coeff-input (list :global global :solution sol-ip)))
	    
	    (when (or (and local-mat (or diffusion-function convection-function))
		      (and local-rhs diffusion-function gamma-function))
	      (let* ((gradients (m* shape-grads Dphi^-1)) ; (n-basis x dim)-matrix
		     (right-gradients (if local-u (m* (transpose local-u) gradients) gradients))
		     (left-gradients (if local-v (m* (transpose local-v) gradients) gradients)))
		;; diffusion
		(when diffusion-function
		  (let* ((diff-tensor (evaluate diffusion-function coeff-input))
			 (fluxes (m* left-gradients diff-tensor))) ; (n-basis x dim)-matrix
		    (when local-mat
		      (gemm! weight fluxes right-gradients 1.0 local-mat :NT))
		    ;; gamma
		    (when (and gamma-function local-rhs)
		      (gemm! weight fluxes (evaluate gamma-function coeff-input)
			     1.0 local-rhs))
		    ))
		;; convection
		(when (and convection-function local-mat)
		  (let* ((velocity-vector (evaluate convection-function coeff-input)))
		    (gemm! (- weight) (m* left-gradients velocity-vector) right-vals 1.0 local-mat :NT)))
	    ))
	    
	    ;; reaction
	    (when (and reaction-function local-mat)
	      (let* ((reaction (evaluate reaction-function coeff-input))
		     (factor (* weight reaction)))
		(when alpha-function
		  (setq factor (* factor (evaluate alpha-function coeff-input))))
		(gemm! factor left-vals right-vals 1.0 local-mat :NT)))
	    
	    ;; source
	    (when (and source-function local-rhs)
	      (let ((source (evaluate source-function coeff-input)))
		(when alpha-function
		  (setq source (scal (evaluate alpha-function coeff-input) source)))
		(when (numberp source) (setq source (make-real-matrix `((,source)))))
		(gemm! weight left-vals source 1.0 local-rhs)
		))
	    ))
    ;; custom fe rhs
    (whereas ((fe-rhs (and local-rhs (getf coeffs 'CDR::FE-RHS))))
      (m+! (funcall fe-rhs (getf fe-geometry :cell) fe) local-rhs))
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
	 (constraints-rhs (make-ansatz-space-vector ansatz-space)))
    (doskel (cell level-skel)
      (when (ecase where
	      (:refined (refined-p cell h-mesh))
	      (:surface (or (not (refined-p cell h-mesh))
			    (member-of-skeleton? cell interface)))
	      (:all t))
	(let* ((coeffs (coefficients-of-cell cell h-mesh problem))
	       (dirichlet-function (getf coeffs 'CDR::DIRICHLET))
	       (cell-key (cell-key cell h-mesh)))
	  (when dirichlet-function
	    (loop with fe = (get-fe fe-class cell)
		  for dof in (fe-dofs fe)
		  for j below (nr-of-inner-dofs fe)
		  for k = (dof-in-vblock-index dof)
		  do
		  (when dirichlet-function
		    ;; The following is only correct for degrees of freedom of
		    ;; Lagrange type.  Perhaps one should use Hermite finite
		    ;; cells only in the interior?
		    (setf (mref (mref constraints-P cell-key cell-key) k k) 1.0)
		    ;;		       (clear-row mat cell k)
		    ;;		       (setf (mref (mref mat cell cell) k k) 1.0)
		    (setf (vref (vref constraints-rhs cell-key) k)
			  (evaluate dirichlet-function
				    (list :local (dof-coord dof)
					  :global (local->global cell (dof-gcoord dof))))))
		  )))))
    (values constraints-P constraints-Q constraints-rhs)))

;;; Testing
(defun cdr-fe-tests ()
  (let* ((order 1) (level 1)
	 (problem (cdr-model-problem 1))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order))
	 (as (make-fe-ansatz-space fedisc problem h-mesh)))
    (with-items (&key matrix rhs)
	(fe-discretize (blackboard :ansatz-space as))
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

