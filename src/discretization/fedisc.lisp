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
;;;; Interface for standard customization (ellsys-fe.lisp)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric discretize-locally (problem coeffs fe qrule fe-geometry
				&key matrix mass-matrix rhs local-u local-v
				fe-parameters)
  (:documentation "Computes a local stiffness matrix and right-hand side.
The algorithm will usually work as follows:

@enumerate
@item Get coefficient functions for the patch of the cell.
@item Compute geometry information for all ips (values and gradients of the shape functions).
@item Loop over integration points ip:
  @enumerate
    @item If necessary, compute input for coefficient functions.
          This input may contain values of finite element function in the
          property list @arg{fe-parameters}.
    @item Evaluate coefficient functions at ips.
    @item Add the contributions for matrix and right-hand side to @arg{local-mat} and @arg{local-rhs}.
  @end enumerate
@end enumerate

If @arg{local-u} and @arg{local-v} are set, then
@arg{local-v}*@arg{local-mat}*@arg{local-u} and
@arg{local-v}*@arg{local-rhs} is computed.  This feature may be used later
on for implementing matrixless computations."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interior assembly
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun assemble-cell (cell ansatz-space &key matrix rhs mass-matrix)
  (let* ((h-mesh (hierarchical-mesh ansatz-space))
         (problem (problem ansatz-space))
         (patch (patch-of-cell cell h-mesh))
         (patch-properties (skel-ref (domain h-mesh) patch)))
    ;; for the moment, we assume that the problem is resolved by the
    ;; domain structure, i.e. that distributional coefficients occur
    ;; only on patches.
    (whereas ((coeffs (filter-applicable-coefficients
                       (coefficients-of-cell cell h-mesh problem)
                       cell patch :constraints nil)))
      (dbg :disc "Coefficients: ~A" coeffs)
      (let* ((fe (get-fe ansatz-space cell))
             (qrule (quadrature-rule fe))
             (geometry (fe-cell-geometry
                        cell (integration-points qrule)
                        :weights (integration-weights qrule)
                        :metric (getf patch-properties 'FL.MESH::METRIC)
                        :volume (getf patch-properties 'FL.MESH::VOLUME)))
             (local-mat (and matrix (make-local-mat ansatz-space cell)))
             (local-mass-mat (and mass-matrix (make-local-mat ansatz-space cell)))
             (local-rhs (and rhs (make-local-vec ansatz-space cell)))
             (fe-paras (loop for obj in (required-fe-functions coeffs)
                             collect obj collect
                                         (get-local-from-global-vec
                                          cell (get-property problem (car obj))))))
        (dbg :disc "FE-parameters: ~A" fe-paras)
        (discretize-locally
         problem coeffs fe qrule geometry
         :matrix local-mat :rhs local-rhs :mass-matrix local-mass-mat
         :fe-parameters fe-paras)
        (dbg :disc "Matrix:~%~A~&Rhs:~%~A" local-mat local-rhs)
        (list :local-mat local-mat :local-rhs local-rhs
              :local-mass-mat local-mass-mat)))))

(defgeneric assemble-interior (ansatz-space where
                               &key level matrix mass-matrix rhs parallel-clustering &allow-other-keys)
  (:documentation
   "Assemble the interior, i.e. ignore constraints arising from boundaries
and hanging nodes.  Discretization is done using the ansatz space
@arg{ansatz-space} on level @arg{level}.  The level argument will usually
be @code{NIL} when performing a global assembly, and be equal to some
number when assembling coarse level matrices for multigrid.  The argument
@arg{where} is a flag indicating where assembly is to be done.  It should
be one of the keywords @code{:surface}, @code{:refined}, @code{:all}.  The
arguments @arg{matrix}, @arg{rhs} should contain vectors/matrices where the
local assembly is accumulated.  Boundary conditions and constraints are not
taken into account within this routine.

In general, this function does most of the assembly work.  Other steps like
handling constraints are intricate, but usually of lower computational
complexity.")
  (:method ((as <ansatz-space>) (where symbol) &rest args &key level parallel-clustering &allow-other-keys)
      (when parallel-clustering (assert (eq where :all)))
    ;; generate list of cells for assembly
    (let* ((cells ())
           (h-mesh (hierarchical-mesh as))
           (skel (if level
                     (cells-on-level h-mesh (if parallel-clustering
                                                (- level parallel-clustering)
                                                level))
                     h-mesh)))
      (doskel (cell skel)
        (when (ecase where
                (:refined (refined-p cell h-mesh))
                (:surface (not (refined-p cell h-mesh)))
                (:all t))
          (push cell cells)))
      (apply #'assemble-interior as cells args)))
  (:method ((as <ansatz-space>) (where list) &key matrix mass-matrix rhs parallel-clustering &allow-other-keys)
      (declare (ignorable parallel-clustering))
    (measure-time-for-block ("~&Ensuring matrix entries needs ~F seconds~%")
      (let ((mesh (mesh matrix)))
        (loop for cell in where
              do
                 (assert (skel-ref mesh cell))
                 (loop
                   with subcells = (subcells cell)
                   for sc1 across subcells
                   for key1 = (cell-key sc1 mesh) do
                     (loop for sc2 across subcells
                           for key2 = (cell-key sc2 mesh)
                           when (entry-allowed-p matrix key1 key2)
                             do (mref matrix key1 key2))))))
    ;; The matrix-transfer-queue contains work accessing the matrix.
    ;; This queue has to be empty, before a new chunk of parallel matrix transfer
    ;; operations can be pushed onto it!
    (let ((chunk-queue (make-instance 'chunk-queue)))
      (flet ((work-on-cell (cell)
               (let ((*use-pool-p* T))
                 ;; try to help other workers with matrix transfer
                 (work-on-queue chunk-queue)
                 ;; then work on assembly
                 (dbg :disc-loop "Working on: ~A" cell)
                 (destructuring-bind (&key local-mat local-rhs local-mass-mat)
                     (assemble-cell cell as :matrix matrix :rhs rhs
                                            :mass-matrix mass-matrix)
                   (flet ((matrix-update (gmat lmat)
                            (when (and gmat lmat)
                              (multiple-value-bind (ops delay-possible-p)
                                  (global-local-matrix-operation
                                   gmat lmat cell cell :global+=local :delay-p t)
                                (cond (delay-possible-p
                                       (chunk-enqueue chunk-queue ops
                                                      (_ (awhen (and *use-pool-p* *local-mat-pool*)
                                                           (put-back-in-pool it lmat))))
                                       ;; and work on matrix transfer
                                       (work-on-queue chunk-queue))
                                      (t 
                                       (with-mutual-exclusion (gmat)
                                         (loop for op in ops do (funcall op)))
                                       (awhen (and *use-pool-p* *local-mat-pool*)
                                         (put-back-in-pool it lmat))))))))
                     ;; accumulate to global matrix and rhs (if not nil)
                     (when (and rhs local-rhs)
                       (increment-global-by-local-vec cell rhs local-rhs))
                     (matrix-update matrix local-mat)
                     (matrix-update mass-matrix local-mass-mat)))
                 )))
        (with-workers (#'work-on-cell :parallel t)  ; !!!
          (dolist (cell where)
            (work-on cell)))
         (when *local-mat-pool*
           (dbg :local-mat-pool "LOCAL-MAT-POOL contains ~D objects after assembly"
                (hash-table-count (slot-value *local-mat-pool* 'fl.parallel::ppes))))
         ))
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Assembly of full problem
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun fe-discretize-linear-problem (blackboard)
  "Finite element discretization for a linear problem or the linearization
of a nonlinear problem."
  (with-items (&key ansatz-space cells (assemble-constraints-p t)
		    matrix mass-matrix rhs solution
		    interior-matrix interior-rhs interior-mass-matrix
		    discretized-problem)
      blackboard
    (let ((problem (problem ansatz-space)))
      (setq interior-matrix (make-ansatz-space-automorphism ansatz-space))
      (setq interior-rhs (make-ansatz-space-vector ansatz-space))
      (when (typep problem '<evp-mixin>)
	(setf interior-mass-matrix (make-ansatz-space-automorphism ansatz-space)))
      (assert (samep (remove nil (list interior-matrix interior-mass-matrix interior-rhs solution))
		      :key #'ansatz-space))
      ;; interior assembly
      (assert (null cells) () "TBI: assembly on cells")
      (assemble-interior ansatz-space :surface
                         :matrix interior-matrix :rhs interior-rhs
                         :mass-matrix interior-mass-matrix)
      (dbg-when :disc
	(dbg :disc "Interior matrix:~%")
	(show interior-matrix)
	(format t "Interior mass matrix:~%")
	(show interior-mass-matrix)
	(format t "Interior rhs:~%")
	(show interior-rhs))
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
	  (setf (get-property result-mat :solution) solution
		(get-property result-rhs :solution) solution)
	  ;; we keep also interior-matrix and interior-rhs which may be of
	  ;; use when assembling for local multigrid.  Note that they will
	  ;; usually share most of their data with result-mat/result-rhs.
	  (setf (get-property result-mat :interior-matrix) interior-matrix
		(get-property result-mat :interior-mass-matrix) interior-mass-matrix
		(get-property result-rhs :interior-rhs)	interior-rhs)
	  (setf matrix result-mat	; might be dropped later on
		mass-matrix interior-mass-matrix
		rhs result-rhs)
	  (assert (and result-mat result-rhs))))))
  ;; return blackboard
  blackboard)

(defun fe-discretize (blackboard)
  "Finite element discretization for an ansatz space provided on the
blackboard."
  (with-items (&key ansatz-space solution matrix mass-matrix rhs discretized-problem)
      blackboard
    (let ((problem (problem ansatz-space)))
      (setf discretized-problem
	    (cond
	      ((typep problem '<evp-mixin>)
	       ;; eigenvalue problem
	       (fe-discretize-linear-problem blackboard)
	       (make-instance '<ls-evp>
			      :stiffness-matrix matrix
			      :mass-matrix mass-matrix
			      :eigenvalues (slot-value problem 'eigenvalues)
			      :multiplicity (multiplicity problem)
			      :solution (or solution (random-ansatz-space-vector ansatz-space))))
	      ((get-property problem 'linear-p)
	       ;; standard linear problem
	       (fe-discretize-linear-problem blackboard)
	       (make-instance '<lse> :matrix matrix :rhs rhs))
	      (t
	       ;; nonlinear problem
	       (make-instance '<nlse>
			      :linearization
			      #'(lambda (sol)
				  (setf (get-property problem :solution) sol)
				  (let ((bb (fe-discretize-linear-problem
					     (blackboard :problem problem
							 :ansatz-space ansatz-space
							 :solution sol))))
				    (with-items (&key matrix rhs) bb
				      (make-instance '<lse> :matrix matrix :rhs rhs))))
			      :solution
			      (or solution (make-ansatz-space-vector ansatz-space))))))))
  blackboard)

(defun discretize-globally (problem h-mesh fe-class)
  "Discretize @arg{problem} on the hierarchical mesh @arg{h-mesh} using
finite elments given by @arg{fe-class}."
  (let ((ansatz-space (make-fe-ansatz-space fe-class problem h-mesh)))
    (with-items (&key matrix rhs)
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
    (return-from discretize (fe-discretize blackboard)))
  (error "You have to provide either an ansatz-space or a mesh in the
blackboard."))

(defun test-fedisc ()
  (let* (;lparallel:*kernel*
         (dim 2) (order 5) (level 6)
         (fe (lagrange-fe order))
         (domain (n-cube-domain dim))
         (problem (cdr-model-problem domain :dirichlet nil))
         (mesh (uniformly-refined-hierarchical-mesh domain level))
         (as (make-fe-ansatz-space fe problem mesh))
         (mat1 (make-ansatz-space-automorphism as))
         (mat2 (make-ansatz-space-automorphism as)))
    (time
     (let (lparallel:*kernel*)
       (assemble-interior as :all :level level :matrix mat1)))
    ;; ensure matrix diagonal of mat2
    (time
     (doskel (cell (cells-on-level mesh level))
       (let ((key (cell-key cell mesh)))
         (when (entry-allowed-p mat2 key key)
           (mref mat2 key key)))))
    (time
     (assemble-interior as :all :level level :matrix mat2 :parallel-clustering 3))
    #+(or)
    (let ((fl.matlisp:*mzerop-threshold* 1e-10))
      (mat-diff mat1 mat2)))
  ;;; (fl.parallel::new-kernel)
  )

