;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  fedisc.lisp
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

(in-package :discretization)

;;; This module contains routines for composing functions for local
;;; assembly to a discretization on a global mesh or for constructing a
;;; multi-level discretization.

;;; It assumes that for the given problem/discretization pair we have
;;; defined methods for the generic functions assemble-interior and
;;; assemble-boundary.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interface for standard customization (e.g. in cdr-fe.lisp)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric discretize-locally (problem coeffs fe qrule fe-geometry
				&key local-mat local-rhs local-sol local-u local-v)
  (:documentation "Computes a local stiffness matrix and right-hand side.
The algorithm will usually work as follows:

- Get coefficient functions for the patch of the cell.

- Loop over integration points ip:
   - Get values and gradients of all shape functions at ip.
   - Compute phi, Dphi, det(Dphi) for the cell transformation phi.
   - Collect necessary information into a coefficient-input structure.
   - Evaluate coefficient functions
   - Compute gradients
   - Compute fluxes.
   - Compute contributions to local-mat/local-rhs for this ip:
     > Multiply fluxes and gradients for diffusion.
     > Multiply source with shape-vals.
     All those contributions are multiplied with quadrature weight and
     volume cell and accumulated into local-mat or local-rhs.

Specialities (not used yet):
- In local-sol is an approximation to the current solution.
- If local-u and/or local-v is provided, the product local-v * local-mat * local-u
  is computed instead of local-mat  
- if local-v is provided, the product local-v * rhs is computed instead of rhs
"))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interior assembly
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun assemble-interior (ansatz-space &key level sol mat rhs (where :surface)
			  &allow-other-keys)
  "Assemble the interior, i.e. ignore constraints arising from boundaries
and hanging nodes.  Parameters are: the ansatz space.  Key parameters are
level, solution, matrix, rhs, and a flag indicating where assembly is to be
done.  This function is used for assembling on each level matrix and rhs,
which are later modified by incorporating boundary conditions and then
combined to yield either a global matrix for use in a direct (sparse)
solver or, alternatively, to assemble levelwise matrices for use within the
multigrid method.

Remark: In this procedure the bulk of the assembly work is done, the
further steps are trickier, but usually of minor complexity."
  (let* ((problem (problem ansatz-space))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (level-skel (if level (cells-on-level h-mesh level) h-mesh))
	 (fe-class (fe-class ansatz-space)))
    (doskel (cell level-skel)
      (let ((patch (patch-of-cell cell h-mesh)))
	;; for the moment, we assume that the problem is resolved by the
	;; domain structure, i.e. that distributional coefficients occur
	;; only on patches.
	(when (= (dimension cell) (dimension patch))
	  (whereas ((coeffs (coefficients-of-cell cell h-mesh problem)))
	    (when (get-properties coeffs (interior-coefficients problem))
	      (when (ecase where
		      (:refined (refined-p cell h-mesh))
		      (:surface (not (refined-p cell h-mesh)))
		      (:all t))
		(let* ((fe (get-fe fe-class cell))
		       (qrule (quadrature-rule fe-class fe))
		       (geometry (fe-cell-geometry cell qrule))
		       (local-mat (and mat (make-local-mat fe)))
		       (local-rhs (and rhs (make-local-vec fe (multiplicity problem))))
		       (local-sol (and sol (get-local-from-global-vec cell fe sol))))
		  (discretize-locally
		   problem coeffs fe qrule geometry
		   :local-mat local-mat :local-rhs local-rhs :local-sol local-sol)
		  ;; accumulate to global matrix and rhs (if not nil)
		  (when rhs (increment-global-by-local-vec cell fe rhs local-rhs))
		  (when mat (increment-global-by-local-mat cell fe mat local-mat)))))))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Assembly of full problem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defun fe-discretize (blackboard)
  "Finite element discretization for an ansatz space provided on the
blackboard."
  (with-items (&key ansatz-space cells (assemble-constraints-p t)
		    matrix rhs solution interior-matrix interior-rhs)
      blackboard
    (setq interior-matrix (or interior-matrix
			      (make-ansatz-space-automorphism ansatz-space)))
    (setq interior-rhs (or interior-rhs
			   (make-ansatz-space-vector ansatz-space)))
    (assert (every #'(lambda (obj) (or (not obj) (eq ansatz-space (ansatz-space obj))))
		   (list interior-matrix interior-rhs solution)))
    ;; interior assembly
    (apply #'assemble-interior ansatz-space :sol solution
	   :mat interior-matrix :rhs interior-rhs
	   (if cells
	       (list :cells cells)
	       '(:where :surface)))
    #+(or)(break)
    (when assemble-constraints-p (assemble-constraints ansatz-space))
    (destructuring-bind (&key constraints-P constraints-Q constraints-r
			      ip-constraints-P ip-constraints-Q ip-constraints-r
			      &allow-other-keys)
	(structure-information ansatz-space)
      #+(or)(break)
      (multiple-value-bind (result-mat result-rhs)
	  (eliminate-constraints interior-matrix interior-rhs
				 ip-constraints-P ip-constraints-Q ip-constraints-r)
	
	;; We put the constraints in the matrix.  An alternative would be
	;; to store the constraints in the ansatz space and to enforce them
	;; after application of the operator
	(m+! constraints-P result-mat)
	(m-! constraints-Q result-mat)
	(m+! constraints-r result-rhs)
	
	#+(or)(break)
	;; we keep also interior-matrix and interior-rhs which may be of
	;; use when assembling for local multigrid.  Note that they will
	;; usually share most of their data with result-mat/result-rhs.
	(setf (getf (discretization-info result-mat) :interior-matrix)
	      interior-matrix)
	(setf (getf (discretization-info result-rhs) :interior-rhs)
	      interior-rhs)
	(setf matrix result-mat
	      rhs result-rhs))))
  ;; return blackboard
  blackboard)

(defun discretize-globally (problem h-mesh fe-class)
  "Old, but effective discretization interface."
  (let ((ansatz-space (make-fe-ansatz-space fe-class problem h-mesh)))
    (destructuring-bind (&key matrix rhs &allow-other-keys)
	(fe-discretize (blackboard :ansatz-space ansatz-space))
      (destructuring-bind (&key ip-constraints-P ip-constraints-Q
				ip-constraints-r &allow-other-keys)
	  (structure-information ansatz-space)
	(values matrix rhs ip-constraints-P ip-constraints-Q ip-constraints-r)))))

(defmethod discretize ((fedisc <fe-discretization>) (problem <problem>) blackboard)
  "General discretization interface for FE."
  (whereas ((as (getbb blackboard :ansatz-space)))
    (assert (and (eq (fe-class as) fedisc)
		 (eq (problem as) problem)))
    (return-from discretize (fe-discretize blackboard)))
  (whereas ((mesh (getbb blackboard :mesh))) 
    (setf (getbb blackboard :ansatz-space)
	  (make-fe-ansatz-space fedisc problem mesh))
    (fe-discretize blackboard))
  (error "You have to provide either an ansatz-space or a mesh in the
blackboard."))
