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

(in-package :fl.discretization)

;;; This module contains routines for composing functions for local
;;; assembly to a discretization on a global mesh or for constructing a
;;; multi-level discretization.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interface for standard customization (e.g. in cdr-fe.lisp)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric discretize-locally (problem coeffs fe qrule fe-geometry
				&key local-mat local-rhs local-sol local-u local-v
				coefficient-parameters)
  (:documentation "Computes a local stiffness matrix and right-hand side.
The algorithm will usually work as follows:

@enumerate
@item Get coefficient functions for the patch of the cell.
@item Compute geometry information for all ips (values and gradients of the shape functions).
@item Loop over integration points ip:
  @enumerate
    @item If necessary, compute input for the coefficient functions.  This input can be another finite element function in the property list @var{coefficient-parameters}.
    @item Evaluate coefficient functions at ips.
    @item Add the contributions for matrix and right-hand side to @var{local-mat} and @var{local-rhs}.
  @end enumerate
@end enumerate

If @var{local-u} and @var{local-v} are set, then
@var{local-v}*@var{local-mat}*@var{local-u} and
@var{local-v}*@var{local-rhs} is computed.  This feature may be used later
on for implementing matrixless computations."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interior assembly
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun assemble-interior (ansatz-space &key level (where :surface)
			  matrix rhs solution &allow-other-keys)
  "Assemble the interior, i.e. ignore constraints arising from boundaries
and hanging nodes.  Parameters are: the ansatz space.  Keyword parameters
are level, solution, matrix, rhs, and a flag indicating where assembly is
to be done.  This function is used for assembling on each level matrix and
rhs, which are later modified by incorporating boundary conditions and then
combined to yield either a global matrix for use in a direct (sparse)
solver or, alternatively, to assemble levelwise matrices for use within the
multigrid method.

In general, this function does most of the assembly work.  Other steps like
handling constraints are tricky, but usually of lower complexity."
  (let* ((problem (problem ansatz-space))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (level-skel (if level (cells-on-level h-mesh level) h-mesh))
	 (fe-class (fe-class ansatz-space)))
    (doskel (cell level-skel)
      (let* ((patch (patch-of-cell cell h-mesh))
	     (patch-properties (skel-ref (domain h-mesh) patch)))
	;; for the moment, we assume that the problem is resolved by the
	;; domain structure, i.e. that distributional coefficients occur
	;; only on patches.
	(when (= (dimension cell) (dimension patch))
	  (whereas ((coeffs (coefficients-of-cell cell h-mesh problem)))
;;;	    (when (get-properties coeffs (interior-coefficients problem))
	    (when (ecase where
		    (:refined (refined-p cell h-mesh))
		    (:surface (not (refined-p cell h-mesh)))
		    (:all t))
	      (let* ((fe (get-fe fe-class cell))
		     (qrule (quadrature-rule fe-class fe))
		     (geometry (fe-cell-geometry
				cell qrule
				:metric (getf patch-properties 'FL.MESH::METRIC)
				:volume (getf patch-properties 'FL.MESH::VOLUME)))
		     (local-mat (and matrix (make-local-mat fe)))
		     (local-rhs (and rhs (make-local-vec fe (multiplicity ansatz-space))))
		     (local-sol (and solution (get-local-from-global-vec cell fe solution))))
		(discretize-locally
		 problem coeffs fe qrule geometry
		 :local-mat local-mat :local-rhs local-rhs :local-sol local-sol)
		;; accumulate to global matrix and rhs (if not nil)
		(when rhs (increment-global-by-local-vec cell fe rhs local-rhs))
		(when matrix (increment-global-by-local-mat cell fe matrix local-mat))))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Assembly of full problem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun fe-discretize-linear-problem (blackboard)
  "Finite element discretization for a linear problem or the linearization
of a nonlinear problem."
  (with-items (&key ansatz-space cells (assemble-constraints-p t)
		    matrix rhs solution interior-matrix interior-rhs
		    discretized-problem)
      blackboard
    (setq interior-matrix (or interior-matrix
			      (make-ansatz-space-automorphism ansatz-space)))
    (setq interior-rhs (or interior-rhs
			   (make-ansatz-space-vector ansatz-space)))
    (assert (every #'(lambda (obj) (or (not obj) (eq ansatz-space (ansatz-space obj))))
		   (list interior-matrix interior-rhs solution)))
    ;; interior assembly
    (assert (null cells) () "TBI: assembly on cells")
    (assemble-interior ansatz-space :solution solution
		       :matrix interior-matrix :rhs interior-rhs
		       :where :surface)
    #+(or)(break)
    (when assemble-constraints-p (assemble-constraints ansatz-space))
    (destructuring-bind (&key constraints-P constraints-Q constraints-r
			      ip-constraints-P ip-constraints-Q ip-constraints-r
			      &allow-other-keys)
	(properties ansatz-space)
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
	
	;; for nonlinear problems the solution vector is important
	(setf (getf (properties result-mat) :solution) solution)
	(setf (getf (properties result-rhs) :solution) solution)
	;; we keep also interior-matrix and interior-rhs which may be of
	;; use when assembling for local multigrid.  Note that they will
	;; usually share most of their data with result-mat/result-rhs.
	(setf (getf (properties result-mat) :interior-matrix)
	      interior-matrix)
	(setf (getf (properties result-rhs) :interior-rhs)
	      interior-rhs)
	(setf matrix result-mat  ; might be dropped later on
	      rhs result-rhs)
	(setf discretized-problem (lse :matrix result-mat :rhs result-rhs)))))
  ;; return blackboard
  blackboard)

(defun fe-discretize (blackboard)
  "Finite element discretization for an ansatz space provided on the
blackboard."
  (with-items (&key ansatz-space solution) blackboard
    (let ((problem (problem ansatz-space)))
      (if (get-property problem 'linear-p)
	  (fe-discretize-linear-problem blackboard)
	  (setf (getbb blackboard :discretized-problem)
		(nlse :linearization
		      #'(lambda (sol)
			  (getbb
			   (fe-discretize-linear-problem
			    (blackboard :problem problem :ansatz-space ansatz-space
					:solution sol))
			   :discretized-problem))))))
    blackboard))

(defun discretize-globally (problem h-mesh fe-class)
  "Discretize @var{problem} on the hierarchical mesh @var{h-mesh} using
finite elments given by @var{fe-class}."
  (let ((ansatz-space (make-fe-ansatz-space fe-class problem h-mesh)))
    (with-items (&key matrix rhs &allow-other-keys)
      (fe-discretize (blackboard :ansatz-space ansatz-space))
      (destructuring-bind (&key ip-constraints-P ip-constraints-Q
				ip-constraints-r &allow-other-keys)
	  (properties ansatz-space)
	(values matrix rhs ip-constraints-P ip-constraints-Q ip-constraints-r)))))

(defmethod discretize ((fedisc <fe-discretization>) (problem <pde-problem>) blackboard)
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
