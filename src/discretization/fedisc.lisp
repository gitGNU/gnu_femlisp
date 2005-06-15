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
				mass-factor stiffness-factor coefficient-parameters)
  (:documentation "Computes a local stiffness matrix and right-hand side.
The algorithm will usually work as follows:

@enumerate
@item Get coefficient functions for the patch of the cell.
@item Compute geometry information for all ips (values and gradients of the shape functions).
@item Loop over integration points ip:
  @enumerate
    @item If necessary, compute input for the coefficient functions.  This input can be another finite element function in the property list @arg{coefficient-parameters}.
    @item Evaluate coefficient functions at ips.
    @item Add the contributions for matrix and right-hand side to @arg{local-mat} and @arg{local-rhs}.
  @end enumerate
@end enumerate

@arg{mass-factor} and @arg{stiffness-factor} are weights for mass and
stiffness matrix which are used for solving eigenvalue and time-dependent
problems.  If @arg{local-u} and @arg{local-v} are set, then
@arg{local-v}*@arg{local-mat}*@arg{local-u} and
@arg{local-v}*@arg{local-rhs} is computed.  This feature may be used later
on for implementing matrixless computations."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interior assembly
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
(defparameter *static-condensation* nil
  "This is experimental and usually switched off, because several things
like error estimators or multigrid do not yet work correctly.")

(defun static-condensation (mat &optional rhs)
  "Does static condensation destructively on mat and rhs, i.e. essentially
eliminate all equations for highest-dimensional cells from the system.  The
system should be completely equivalent to the original one."
  (declare (optimize debug safety))
  (let ((dim (dimension (mesh mat))))
    (for-each-row-key
     #'(lambda (i)
	 (when (= (dimension (representative i)) dim)
	   (let ((Aii-inv (m/ (mref mat i i)))
		 (f_i (when rhs (vref rhs i))))
	     (for-each-key-in-col
	      #'(lambda (l)
		  (unless (eq l i)
		    (let ((B (m* (mref mat l i) Aii-inv)))
		      (when rhs (gemm! -1.0 B  f_i 1.0 (vref rhs l)))
		      (for-each-key-in-row
		       #'(lambda (j)
			   (unless (eq j i)
			     (gemm! -1.0 B (mref mat i j) 1.0 (mref mat l j))))
		       mat i))
		    (remove-entry mat l i)))
	      mat i))))
     mat)))

(defun assemble-interior (ansatz-space &key level (where :surface)
			  matrix rhs solution ; left right
			  (mass-factor 0.0) (stiffness-factor 1.0)
			  &allow-other-keys)
  "Assemble the interior, i.e. ignore constraints arising from boundaries
and hanging nodes.  Discretization is done using the ansatz space
@arg{ansatz-space} on level @arg{level}.  The level argument will usually
be @code{NIL} when performing a global assembly, and be equal to some
number when assembling coarse level matrices for multigrid.  The argument
@arg{where} is a flag indicating where assembly is to be done.  It should
be one of the keywords @code{:surface}, @code{:refined}, @code{:all}.  The
arguments @arg{solution}, @arg{matrix}, @arg{rhs} should contain
vectors/matrices where the local assembly is accumulated.  The numbers
@arg{mass-factor} and @arg{stiffness-factor} determine weights for mass and
stiffness matrix which is used when solving time-dependent and eigenvalue
problems.  Boundary conditions and constraints are not taken into account
within this routine.

In general, this function does most of the assembly work.  Other steps like
handling constraints are intricate, but usually of lower complexity."
  (dbg :disc "Mass-factor=~A Stiffness-factor=~A" mass-factor stiffness-factor)
  (let* ((problem (problem ansatz-space))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (level-skel (if level (cells-on-level h-mesh level) h-mesh))
	 (fe-class (fe-class ansatz-space)))
    (doskel (cell level-skel)
      (when (ecase where
	      (:refined (refined-p cell h-mesh))
	      (:surface (not (refined-p cell h-mesh)))
	      (:all t))
	(let* ((patch (patch-of-cell cell h-mesh))
	       (patch-properties (skel-ref (domain h-mesh) patch)))
	  ;; for the moment, we assume that the problem is resolved by the
	  ;; domain structure, i.e. that distributional coefficients occur
	  ;; only on patches.
	  (whereas ((coeffs (filter-applicable-coefficients
			     (coefficients-of-cell cell h-mesh problem)
			     cell patch :constraints nil)))
	    (let* ((fe (get-fe fe-class cell))
		   (qrule (quadrature-rule fe))
		   (geometry (fe-cell-geometry
			      cell qrule
			      :metric (getf patch-properties 'FL.MESH::METRIC)
			      :volume (getf patch-properties 'FL.MESH::VOLUME)))
		   (local-mat (and matrix (make-local-mat fe)))
		   (local-rhs (and rhs (make-local-vec fe (multiplicity ansatz-space))))
		   (local-sol (and solution (get-local-from-global-vec cell fe solution))))
	      (discretize-locally
	       problem coeffs fe qrule geometry
	       :local-mat local-mat :local-rhs local-rhs :local-sol local-sol
	       :mass-factor mass-factor :stiffness-factor stiffness-factor)
	      ;; accumulate to global matrix and rhs (if not nil)
	      (when rhs (increment-global-by-local-vec cell fe rhs local-rhs))
	      (when matrix (increment-global-by-local-mat cell fe matrix local-mat))))))))
  ;;  experimental: static condensation
  (when *static-condensation*
    (static-condensation matrix rhs))
  )

(defun compute-interior-level-matrix (interior-mat sol level)
  "This function is needed for the multilevel decomposition of geometric
multigrid."
  (let* ((ansatz-space (ansatz-space interior-mat))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (top-level (top-level h-mesh))
	 (mat (extract-level interior-mat level)))
    ;; extend the surface matrix on this level by the refined region
    (when (< level top-level)
      (assemble-interior ansatz-space :matrix mat :solution sol :level level :where :refined))
    ;; extend it by the hanging-node region
    (loop for level from (1- level) downto 0
	  for constraints =
	  (nth-value 1 (hanging-node-constraints ansatz-space :level level :ip-type t))
	  for constraints-p =
	  (eliminate-hanging-node-constraints-from-matrix mat constraints)
	  while constraints-p do
	  (add-local-part! mat interior-mat (column-table constraints)
			   :directions '(:right)))
    ;; experimental: static condensation
    (when *static-condensation*
      (static-condensation mat))
    ;; and return the result
    mat))


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
    (let ((problem (problem ansatz-space))
	  (mass-factor 0.0)
	  (stiffness-factor 1.0))
      (typecase problem
	(<evp-mixin>
	 (setf mass-factor (- (unbox (slot-value problem 'lambda)))
	       stiffness-factor (unbox (slot-value problem 'mu)))))
      (setq interior-matrix (make-ansatz-space-automorphism ansatz-space))
      (setq interior-rhs (make-ansatz-space-vector ansatz-space))
      (assert (every #'(lambda (obj) (or (not obj) (eq ansatz-space (ansatz-space obj))))
		     (list interior-matrix interior-rhs solution)))
      ;; interior assembly
      (assert (null cells) () "TBI: assembly on cells")
      (assemble-interior ansatz-space :solution solution
			 :matrix interior-matrix :rhs interior-rhs
			 :where :surface
			 :mass-factor (force mass-factor)
			 :stiffness-factor (force stiffness-factor))
      (when assemble-constraints-p (assemble-constraints ansatz-space))
      (destructuring-bind (&key constraints-P constraints-Q constraints-r
				ip-constraints-P ip-constraints-Q ip-constraints-r
				&allow-other-keys)
	  (properties ansatz-space)
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
	  (setf matrix result-mat	; might be dropped later on
		rhs result-rhs)
	  (assert (and result-mat result-rhs))
	  (setf discretized-problem (lse :matrix result-mat :rhs result-rhs))))))
  ;; return blackboard
  blackboard)

(defun fe-discretize (blackboard)
  "Finite element discretization for an ansatz space provided on the
blackboard."
  (with-items (&key ansatz-space solution)
      blackboard
    (let ((problem (problem ansatz-space)))
      (if (get-property problem 'linear-p)
	  (fe-discretize-linear-problem blackboard)
	  (setf (getbb blackboard :discretized-problem)
		(let ((linearization
		       #'(lambda (sol)
			   (getbb (fe-discretize-linear-problem
				   (blackboard :problem problem :ansatz-space ansatz-space :solution sol))
				  :discretized-problem))))
		  (if (typep problem '<evp-mixin>)
		      (make-instance '<evp> :linearization linearization
				     :lambda (slot-value problem 'lambda)
				     :mu (slot-value problem 'mu)
				     :initial-guess
				     (or solution
					 (random-ansatz-space-vector ansatz-space)))
		      (make-instance '<nlse> :linearization linearization
				     :initial-guess
				     (or solution
					 (make-ansatz-space-vector ansatz-space)))))))))
  blackboard)

(defun discretize-globally (problem h-mesh fe-class)
  "Discretize @arg{problem} on the hierarchical mesh @arg{h-mesh} using
finite elments given by @arg{fe-class}."
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
