;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; geomg.lisp
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

(in-package :geomg)

;;; This file provides the geometric multigrid iteration which is a
;;; multigrid iteration working with hierarchical ansatz-space vectors
;;; coming from discretizations on refined grids.  At the moment, our
;;; version works only on uniformly refined grids.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Geometric multigrid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <geometric-mg> (<mg-iteration>)
  ((solution :accessor solution :initarg :solution :type <ansatz-space-vector>))
  (:documentation "The geometric multigrid iteration is a multigrid
iteration where the hierarchy of problems is obtained by discretizing on a
sequence of refined meshes.  It is an abstract class and should be merged
with either <correction-scheme> or <fas>."))

(defun extend-horizontally (mat surface-mat)
  "Adds one layer from surface-mat to mat."
  (for-each-col-key
   #'(lambda (key)
       (unless (matrix-row mat key)
	 (for-each-key-and-entry-in-row
	  #'(lambda (ck entry)
	      (setf (mat-ref mat key ck) entry))
	  surface-mat key)))
   mat))

(defun extend-level-matrix (initial-mat surface-mat &key (layers 0))
  "Extends initial-mat by several layers from surface-mat."
  (destructuring-bind (&key hanging-P hanging-Q hanging-r
			    &allow-other-keys)
      (structure-information (ansatz-space initial-mat))
    (loop with mat = (eliminate-constraints
		      initial-mat nil hanging-P hanging-Q hanging-r
		      :assemble-locally t :include-constraints t)
	  repeat layers do
	  (extend-horizontally mat surface-mat)
	  (setq mat (eliminate-constraints
		     mat nil hanging-P hanging-Q hanging-r
		     :assemble-locally t :include-constraints t))
	  finally (return mat))))

(defmethod multilevel-decomposition ((mgit <geometric-mg>)
				     (mat <ansatz-space-automorphism>))
  "Assemble all levels below the top-level.  The top level should have been
already assembled.  Works only for uniformly refined meshes."
  (let* ((ansatz-space (ansatz-space mat))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (top-level (top-level h-mesh))
	 (a-vec (make-array (nr-of-levels h-mesh)))
	 (i-vec (make-array (nr-of-levels h-mesh)))
	 (interior-mat (getf (discretization-info mat) :interior-matrix)))
    
    ;; set the matrix vector
    (multiple-value-bind (essential-P essential-Q essential-r)
	(discretization::compute-essential-boundary-constraints
	 ansatz-space :where :all)
      (loop for level from 0 upto top-level
	    for l-mat = (discretization::compute-interior-level-matrix
			 interior-mat level) do
	    (let ((eliminated-mat
		   (eliminate-constraints
		    l-mat nil essential-P essential-Q essential-r
		    :include-constraints t)))
	      (setf (aref a-vec level) eliminated-mat))))
    ;; set the interpolation vector
    (loop for level below top-level
	  for imat = (constrained-interpolation-matrix
		      ansatz-space :level level :where :refined)
	  do
	  (extend-by-identity imat (row-table (aref a-vec level))
			      :ignore (column-table imat))
	  (setf (aref i-vec level) imat))
    ;; return result
    (blackboard :a-vec a-vec :i-vec i-vec)))

#+(or)
(defmethod multilevel-decomposition ((mgit <geometric-mg>)
				     (mat <ansatz-space-automorphism>))
  "Assemble all levels below the top-level.  The top level should have been
already assembled.  Works only for uniformly refined meshes."
  (let* ((ansatz-space (ansatz-space mat))
	 (h-mesh (hierarchical-mesh ansatz-space))
	 (top-level (top-level h-mesh))
	 (a-vec (make-array (nr-of-levels h-mesh)))
	 (i-vec (coerce
		 (loop for level below (top-level (mesh ansatz-space))
		       collect (constrained-interpolation-matrix
				ansatz-space :level level :where :refined))
		 'vector)))
    ;; set the matrix vector
    (setf (aref a-vec top-level) mat)
    (loop for level below top-level do
	  (let ((l-mat (make-ansatz-space-automorphism ansatz-space)))
	    (assemble-interior ansatz-space :mat l-mat :level level :where :refined)
	    (multiple-value-bind (constraints-P constraints-Q constraints-r)
		(discretization::compute-essential-boundary-constraints
		 ansatz-space :level level :where :refined)
	      (let ((eliminated-mat (eliminate-constraints
				     l-mat nil constraints-P constraints-Q constraints-r)))
		(x+=y eliminated-mat constraints-P)
		(x-=y eliminated-mat constraints-Q)
		(setf (aref a-vec level) eliminated-mat)))))
    ;; return result
    (make-instance '<mg-data> :a-vec a-vec :i-vec i-vec)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; correction scheme
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <geometric-cs> (<correction-scheme> <geometric-mg>)
  ()
  (:documentation "Geometric multigrid of correction scheme type."))

(defun geometric-cs (&rest key-args)
  "Constructor of a geometric multigrid iteration of correction scheme
type."
  (apply #'make-instance '<geometric-cs> key-args))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; local multigrid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <local-mg> (<geometric-cs>)
  ()
  (:documentation "Local geometric multigrid of correction scheme type."))

(defun local-mg (&rest key-args)
  "Constructor of a geometric multigrid iteration of correction scheme
type."
  (apply #'make-instance '<local-mg> key-args))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; FAS (full approximation scheme)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <geometric-fas> (<fas> <geometric-mg>)
  ()
  (:documentation "Brandt's FAS scheme approximates the unknowns on every
level instead of using corrections.  This requires slightly more work, but
is better suited for handling nonlinear problems and local refinements."))

(defun fas (&rest key-args)
  "Constructor of a geometric multigrid iteration of FAS type."
  (apply #'make-instance '<geometric-fas> key-args))

(defmethod multilevel-decomposition ((mgit <geometric-fas>) (mat <ansatz-space-automorphism>))
  "Add the vector of FAS restrictions to the output of this method for
<geometric-mg>."
  (let ((mg-data (call-next-method)))
    (setf (getbb mg-data :fas-r-vec)
	  (let ((ansatz-space (ansatz-space mat)))
	    (coerce
	     (loop for level from 0 below (top-level (mesh ansatz-space))
		   collect (projection-matrix ansatz-space :level level))
	     'vector)))
    mg-data))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <s1-reduction>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The following is a utility which can make scalar higher order problems
;;; treatable by AMG.

(defclass <s1-reduction> (<selection-amg>)
  ()
  (:documentation "This class is useful for reducing a higher-order FE
discretization to a first-order FE discretization.  This can afterwards be
treated by ordinary AMG steps.  Even if it has the structure of a
<selection-amg>, it is far from being a pure algebraic multigrid."))

#-(or)
(defmethod prolongation ((amg <s1-reduction>) (mat <ansatz-space-automorphism>))
  "Transfer to S1 finite elements."
  (let* ((as (ansatz-space mat))
	 (as-1 (make-fe-ansatz-space
		(lagrange-fe 1) (problem as) (mesh as))))
    (assemble-constraints as-1)
    (transfer-matrix as-1 as :no-slaves t)))

#+(or)  ; old: will be dropped soon
(defmethod prolongation ((amg <s1-reduction>) (mat <ansatz-space-automorphism>))
  "This prolongation is somewhat tricky, because it includes also the case
of hanging degrees of freedom and does not interpolate from slaves."
  (declare (optimize (debug 3)))
  (let* ((ansatz-space (ansatz-space mat))
	 (constraints-P (getf (discretization::structure-information ansatz-space)
			      :ip-constraints-P))
	 (constraints-Q (getf (discretization::structure-information ansatz-space)
			      :ip-constraints-Q))
	 (mesh (mesh ansatz-space))
	 (lagrange-1 (lagrange-fe 1))
	 ;;(domain-as (make-fe-ansatz-space lagrange-1 (problem mat) mesh))
	 (prol (make-ansatz-space-automorphism ansatz-space)))
    ;;(assert (symmetric-p mat :threshold 1.0e-7 :output t))  ;;!!
    ;; first pass: represent non-vertex-dofs by vertex-dofs
    (for-each-row-key
     #'(lambda (key)
	 (assert (not (slave-or-dirichlet-dof-p key mat))) ; should be filtered out
	 (let ((cell (representative key)))
	   (let* ((image-fe (get-fe (fe-class ansatz-space) cell))
		  (domain-fe (get-fe lagrange-1 cell))
		  (local-prol (local-transfer-matrix domain-fe image-fe))
		  (subcells (subcells cell))
		  (nrows (nr-of-inner-dofs image-fe)))
	     (loop for subcell across subcells
		   for subcell-ndofs across (subcell-ndofs domain-fe)
		   and l = 0 then (+ l subcell-ndofs)
		   unless (zerop subcell-ndofs) do
		   (let ((subcell-key (cell-key subcell mesh)))
		     (assert (vertex? (representative subcell-key)))
		     (setf (mat-ref prol key subcell-key)
			   (matrix-slice
			    local-prol
			    :from-row 0 :from-col l
			    :nrows nrows :ncols subcell-ndofs)))))))
     mat)
    ;; next pass: remove Dirichlet dofs
    (loop for key in (col-keys prol)
	  when (c-dirichlet-dof-p key constraints-P constraints-Q)
	  do #+debug(format t "~&remove Dirichlet ~A~%" key)
	  (remove-column prol key))
    ;; next pass: eliminate hanging nodes and non-vertices, maybe recursively
    (let ((problem-keys (remove-if-not (rcurry #'c-slave-dof-p constraints-P constraints-Q)
				       (col-keys prol))))
      (loop
       while problem-keys do
       (loop with problem-table = (make-hash-table)
	     for j in problem-keys
	     for slave-p = (c-slave-dof-p j constraints-P constraints-Q)
	     do
	     #+debug (format t "Handling problem key ~A slave=~A~%" j slave-p)
	     (for-each-key-and-entry-in-col
	      #'(lambda (i P_ij)
		  ;; eliminate P_ij (replace j in terms of its representation)
		  #+debug (format t "   column ~A~%" i)
		  (for-each-key-and-entry-in-row
		   #'(lambda (k X_jk)
		       (assert (not (eq j k)))
		       (unless (eq j k)
			 (gemm! 1.0d0 P_ij X_jk 1.0d0 (mat-ref prol i k))
			 #+debug (format t "     substituting ~A ~A~%" k X_jk)
			 (when (or (c-slave-dof-p k constraints-P constraints-Q)
				   (not (vertex? (representative k))))
			   (setf (gethash k problem-table) t))))
		   (if slave-p constraints-Q prol) j))
	      prol j)
	     (remove-column prol j)
	     finally (setq problem-keys (hash-table-keys problem-table)))))
    #+(or)
    (mat-diff prol
	      (let* ((as (ansatz-space mat))
		     (as-1 (make-fe-ansatz-space
			    (lagrange-fe 1) (problem as) (mesh as))))
		(assemble-constraints as-1)
		(transfer-matrix as-1 as)))
    prol))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; our default amg solver for higher-order equations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <s1-coarse-grid-iterator> (<linear-iteration>)
  ()
  (:documentation "Calls the coarse-grid solver of *standard-stueben*
directly, if the matrix was not reduced to S1 which may happen if there are
only Dirichlet vertices."))

(defmethod make-iterator ((s1-cgit <s1-coarse-grid-iterator>) (A <sparse-matrix>))
  (let ((all-vertices-p t))
    (for-each-row-key #'(lambda (key)
			  (unless (vertex? (representative key))
			    (setq all-vertices-p nil)))
		      A)
    (assert (symmetric-p A :threshold 1.0e-5))
    (if all-vertices-p
	(make-iterator (make-instance '<stueben>) A)
	(make-iterator (coarse-grid-iteration *standard-stueben*) A))))

(defun s1-reduction-amg-solver (order &key reduction output (maxsteps 100))
  "This is an AMG solver which works also for Lagrange fe of order p by
reducing them to P^1 first."
  (declare (ignore order))
  (let ((smoother (geometric-ssc)))
    (make-instance
     '<linear-solver>
     :iteration
     (make-instance '<s1-reduction> :max-depth 2 :pre-steps 1 :pre-smooth smoother
		    :post-steps 1 :post-smooth smoother
		    :coarse-grid-iteration #+(or)*lu-iteration*
		    #-(or)(make-instance '<s1-coarse-grid-iterator>)
		    :output output)
     :success-if `(< :reduction ,(or reduction 1.0e-2))
     :failure-if `(> :step ,maxsteps)
     :output output)))


;;; Testing: (test-geomg)

(defun test-geomg ()
  (geometric-cs :fmg t)
  (fas :fmg t :coarse-grid-iteration *lu-iteration*)
  (s1-reduction-amg-solver 2)
  (class-name (class-of (make-instance '<s1-reduction>)))
  )

(tests::adjoin-femlisp-test 'test-geomg)
