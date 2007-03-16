;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; fe-transfer.lisp
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

;;;; This file provides routines for computing local transfer routines
;;;; between fes for a cell and its refinement as well as between fes on
;;;; the same cell.

(defparameter *ipt-tolerance* 1.0e-10
  "Entries in interpolation, projection, and transfer matrices below this
threshold are dropped.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun children-offsets (vecfe type rule)
  "Computes the offsets for each component in the children's data block."
  (let* ((nr-comps (nr-of-components vecfe))
	 (refcell (reference-cell vecfe))
	 (children (ecase type
		     (:inner (inner-refcell-children refcell rule))
		     (:all (refcell-children refcell rule))))
	 (result (make-array nr-comps :initial-element nil)))
    (dotimes (comp nr-comps result)
      (setf (aref result comp)
	    (map 'vector
		 #'(lambda (child)
		     (let ((child-fe (get-fe (discretization vecfe) child)))
		       (aref (aref (subcell-offsets child-fe) comp) 0)))
		 children)))))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interpolation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric compute-local-imatrix (rule fe child-disc &key &allow-other-keys)
  (:documentation "Returns a local interpolation matrix as a sparse matrix
with standard-matrix entries."))

(defmethod compute-local-imatrix
    (rule (fe <scalar-fe>) (child-disc <scalar-fe-discretization>)
     &key (shape-distortion #'identity) (type :local))
  "This function computes a local interpolation matrix for the coarse
finite element @arg{fe}, the refinement rule @arg{rule} and the
discretization @arg{child-disc} of the children.  For each refined interior
degree of freedom, it evaluates the corresponding degree of freedom on the
basis polynomials of the father cell which yields the entry in the local
interpolation matrix.  @arg{shape-distorion} allows to modify the shape
which can be used for implementing problem-dependent finite elements."
  (let* ((refcell (reference-cell fe))
	 (children
	  (case type
	    (:local (inner-refcell-children refcell rule))
	    (t (refcell-children refcell rule))))
	 (subcells (subcells refcell))
	 (subcell-ndofs (subcell-ndofs fe))
	 (imat (make-instance '<sparse-tensor> :rank 2)))
    (loop for child across children and i from 0
	  for child-fe = (get-fe child-disc child)
	  when (plusp (nr-of-inner-dofs child-fe)) do
	  (loop for j below (length subcells)
	     when (plusp (aref subcell-ndofs j)) do
	       (let ((mblock (make-real-matrix (nr-of-inner-dofs child-fe)
					       (aref subcell-ndofs j))))
		 (do-dof (child-dof child-fe)
		   (when (interior-dof? child-dof)
		     (do-dof ((dof shape) fe :type :dof-and-shape)
		       (when (= (dof-subcell-index dof) j)
			 (setf (mref mblock (dof-in-vblock-index child-dof)
				     (dof-in-vblock-index dof))
			       (evaluate child-dof
					 (compose-2 (funcall shape-distortion shape)
						    (cell-mapping child))))))))
		 (unless (mzerop mblock *ipt-tolerance*)
		   (setf (tensor-ref imat i j) mblock)))))
    imat))

(defmethod compute-local-imatrix
    (rule (vecfe <vector-fe>) (child-disc <vector-fe-discretization>)
     &key (type :local) &allow-other-keys)
  "Computes an interpolation matrix for the refinement rule @arg{rule} and
the vector finite element @arg{vecfe}.  The algorithm evaluates the nodal
functionals of the children on the parent shape functions."
  (assert (eq type :local))
  (with-slots (refcell components properties)
    vecfe
    (let* ((children (inner-refcell-children refcell rule))
	   (subcells (subcells refcell))
	   (subcell-ndofs (subcell-ndofs vecfe))
	   (subcell-offsets (subcell-offsets vecfe))
	   (children-offsets (children-offsets vecfe :inner rule))
	   (fe-imats (vector-map (curry #'compute-local-imatrix rule)
				 components (components child-disc)))
	   (vecfe-imat (make-instance '<sparse-tensor> :rank 2)))
      (loop for child across children and i from 0
	    for child-fe = (get-fe child-disc child) do
	    (loop for j below (length subcells)
	       when (some #'(lambda (imat) (in-pattern-p imat i j))
			  fe-imats) do
		 (loop with mblock = (make-real-matrix (nr-of-inner-dofs child-fe)
							(aref subcell-ndofs j))
		    for fe-imat across fe-imats and comp from 0
		    for row-off = (aref (aref children-offsets comp) i)
		    and col-off = (aref (aref subcell-offsets comp) j)
		    when (in-pattern-p fe-imat i j) do
		      (minject (tensor-ref fe-imat i j) mblock row-off col-off)
		    finally (setf (tensor-ref vecfe-imat i j) mblock))))
      vecfe-imat)))

(with-memoization (:id 'local-imatrix)
  (defun local-imatrix (rule fe &optional (type :local))
    "Memoized call of compute-local-imatrix."
    (assert (discretization fe))
    (assert (eq (reference-cell rule) (reference-cell fe)))
    (dbg :fe "Generating imatrix for rule ~A and fe ~A" rule fe)
    (memoizing-let ((rule rule) (disc (discretization fe)) (type type))
      (compute-local-imatrix rule fe disc :type type))))

(defvar *local-interpolation-matrix* nil
  "If non-NIL, it should be a function of the arguments cell, mesh, and
finite element class.  This function is called for computing an
interpolation matrix.")

(defun local-interpolation-matrix (cell ansatz-space type)
  "Returns a local interpolation matrix for interpolating the given
@arg{fe-class} from @arg{cell} to its children in @arg{h-mesh}.  If
@arg{type} is :local, interpolation extends only to the interior children."
  (aif *local-interpolation-matrix*
       (funcall it cell ansatz-space type)
       (local-imatrix
	(refinement-rule cell (hierarchical-mesh ansatz-space))
	(get-fe ansatz-space cell)
	type)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Projection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric compute-local-pmatrix (rule fe child-disc)
  (:documentation "Computes a local projection matrix as a sparse vector
with standard-matrix entries.  Note that projection is not as canonic as
interpolation.  We implement here the injection of dof functionals,
i.e. each degree of freedom is evaluated for the refined fe function."))

(defmethod compute-local-pmatrix (rule (fe <scalar-fe>) child-disc)
  "This function computes a local projection matrix for Lagrangian finite
cells and maybe others.  For each interior degree of freedom, it finds the
child cell where the support of the node functional is contained and
evaluates it on the basis polynomials of the child which yields the entry
in the local projection matrix."
  (let ((nr-parent-inner-dofs (nr-of-inner-dofs fe)))
    (when (plusp nr-parent-inner-dofs)
      (let* ((refcell (reference-cell fe))
	     (all-children (refcell-children refcell rule))
	     (pmat (make-instance '<sparse-tensor> :rank 1))
	     (parent-dof-index 0))
	(do-dof (dof fe)
	  (when (interior-dof? dof)
	    (let ((child (find-if #'(lambda (child)
				      (and (= (dimension child) (dimension refcell))
					   (inside-cell? child (dof-coord dof))))
				  all-children)))
	      (assert child () "No child for dof ~A found" dof)
	      (let ((child-fe (get-fe child-disc child))
		    (child-subcells (subcells child)))
		(do-dof ((child-dof child-shape) child-fe :type :dof-and-shape)
		  (let ((entry (evaluate dof (compose-2 child-shape (curry #'global->local child)))))
		    (when (> (abs entry) *ipt-tolerance*)
		      (let* ((subchild-index (dof-subcell-index child-dof))
			     (subcell (aref child-subcells subchild-index))
			     (subchild-pos (position subcell all-children))
			     (subchild-ndofs (aref (subcell-ndofs child-fe) subchild-index)))
			(assert subchild-pos)
			(unless (in-pattern-p pmat subchild-pos)
			  (setf (tensor-ref pmat subchild-pos)
				(make-real-matrix nr-parent-inner-dofs subchild-ndofs)))
			(setf (mref (tensor-ref pmat subchild-pos)
				    parent-dof-index (dof-in-vblock-index child-dof))
			      entry))))))))
	  (incf parent-dof-index))
	;; return pmat
	pmat))))

(defmethod compute-local-pmatrix (rule (vecfe <vector-fe>) child-disc)
  "Computes a local projection matrix for vector finite elements."
  (declare (optimize (debug 3)))
  (with-slots (refcell discretization components properties)
      vecfe
    (let ((nr-parent-inner-dofs (nr-of-inner-dofs vecfe)))
      (when (plusp nr-parent-inner-dofs)
	(let* ((all-children (refcell-children refcell rule))
	       (children-offsets (children-offsets vecfe :all rule))
	       (subcell-offsets (subcell-offsets vecfe))
	       (fe-pmats (vector-map (curry #'compute-local-pmatrix rule)
				     components (components child-disc)))
	       (vecfe-pmat (make-instance '<sparse-tensor> :rank 1)))
	  (loop for child across all-children and i from 0
	     for child-fe = (get-fe child-disc child)
	     when (some #'(lambda (pmat) (and pmat (in-pattern-p pmat i))) fe-pmats) do
	     (loop with mblock = (make-real-matrix nr-parent-inner-dofs
						   (nr-of-inner-dofs child-fe))
		for pmat across fe-pmats and comp from 0
		for child-off = (aref (aref children-offsets comp) i)
		for parent-off = (aref (aref subcell-offsets comp) 0)
		when (and pmat (in-pattern-p pmat i)) do
		  (minject (tensor-ref pmat i) mblock parent-off child-off)
		finally (setf (tensor-ref vecfe-pmat i) mblock)))
	  vecfe-pmat)))))

(with-memoization (:id 'local-pmatrix)
  (defun local-pmatrix (rule fe)
    "Memoized call of compute-local-pmatrix."
    (assert (eq (reference-cell rule) (reference-cell fe)))
    (memoizing-let ((rule rule) (fe fe))
      (dbg :fe "Generating projection matrix for rule ~A and fe ~A" rule fe)
      (compute-local-pmatrix rule fe (discretization fe)))))

(defvar *local-projection-matrix* nil
  "If non-NIL, it should be a function of the arguments cell, mesh, and
finite element class.  This function is called for computing a projection
matrix.")

(defun local-projection-matrix (cell ansatz-space)
  "Returns a local projection matrix for projecting the given
@arg{fe-class} from the children of @arg{cell} in @arg{h-mesh} to
@arg{cell}."
  (aif *local-projection-matrix*
       (funcall it cell ansatz-space)
       (local-pmatrix (refinement-rule cell (hierarchical-mesh ansatz-space))
		      (get-fe ansatz-space cell))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Transfer between different fe-spaces
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(with-memoization ()
  (defun local-transfer-matrix (fe-from fe-to)
    "Computes a local transfer matrix between different FE spaces."
    (assert (= 1 (nr-of-components fe-from) (nr-of-components fe-to)))
    (memoizing-let ((fe-from fe-from) (fe-to fe-to))
      (assert (eq (reference-cell fe-from) (reference-cell fe-to)))
      (let ((local-mat (make-real-matrix (nr-of-dofs fe-to) (nr-of-dofs fe-from))))
	(loop+ (i (dof (fe-dofs fe-to)))
	  do (loop+ (j (phi (fe-basis fe-from)))
	       do (setf (mref local-mat i j) (evaluate dof phi))))
	local-mat))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing (test-fetransfer)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-fetransfer ()
  (flet ((refrule (dim)
	   (get-refinement-rule (n-cube dim) :regular))
	 (fe (order &optional comp (dim 1))
	   (get-fe (lagrange-fe order :nr-comps comp)
		   (n-cube dim))))
    (fe 1)
    (local-imatrix (refrule 1) (fe 1))
    (let ((local-imat (local-imatrix (refrule 1) (fe 1))))
      (show local-imat)
      (tensor-ref local-imat 0 1))
    (show (local-imatrix (refrule 1) (fe 3)))
    (show (local-pmatrix (refrule 1) (fe 3)))
    (show (local-imatrix (refrule 0) (fe 2 2 0)))
    (show (local-imatrix (refrule 1) (fe 2 2)))
    (assert
     (= 28 (total-entries (local-imatrix (refrule 1) (fe 2 2)))))
    (show (local-pmatrix (refrule 1) (fe 4)))
    (show (local-pmatrix (refrule 0) (fe 1 1 0)))
    (show (local-pmatrix (refrule 0) (fe 1 2 0)))
    (show (local-pmatrix (refrule 1) (fe 4 2)))
    (let ((count 0))
      (dotensor (entry (local-imatrix (refrule 1) (fe 1)) :depth 2)
	(when entry (incf count)))
      (assert (= count 2)))
    (describe (tensor-ref (local-imatrix (refrule 1) (fe 1)) 0))
    (let* ((fedisc (lagrange-fe 1))
	   (fe (get-fe fedisc (n-cube 1))))
      (compute-local-imatrix
       :regular fe fedisc
       :shape-distortion #'(lambda (poly) (poly-expt poly 2))))
    ))

;;; (test-fetransfer)
(fl.tests:adjoin-test 'test-fetransfer)
