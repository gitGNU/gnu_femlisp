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

(defun children-offsets (vecfe type)
  "Computes the offsets for each component in the children's data block."
  (let* ((nr-comps (nr-of-components vecfe))
	 (refcell (reference-cell vecfe))
	 (children (ecase type
		     (:inner (inner-refcell-children refcell))
		     (:all (refcell-children refcell))))
	 (result (make-array nr-comps)))
    (dotimes (comp nr-comps result)
      (setf (aref result comp)
	    (let ((child-offsets (make-array (length children))))
	      (dotimes (j (length children) child-offsets)
		(let* ((child (aref children j))
		       (child-fe (get-fe (discretization vecfe) child)))
		  (setf (aref child-offsets j)
			(aref (aref (subcell-offsets child-fe) comp) 0)))))))))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interpolation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric compute-local-imatrix (fe)
  (:documentation "Returns a local interpolation matrix as a sparse matrix
with standard-matrix entries."))
  
(defmethod compute-local-imatrix ((fe <fe>))
  "This function computes a local interpolation matrix for scalar
Lagrangian finite elements and maybe other types as well.  For each refined
interior degree of freedom, it evaluates the corresponding degree of
freedom on the basis polynomials of the father cell which yields the entry
in the local interpolation matrix."
  (let* ((refcell (reference-cell fe))
	 (fe-disc (discretization fe))
	 (children (inner-refcell-children refcell))
	 (subcells (subcells refcell))
	 (subcell-ndofs (subcell-ndofs fe))
	 (imat (make-instance '<sparse-tensor> :rank 2)))
    (loop for child across children and i from 0
	  for child-fe = (get-fe fe-disc child)
	  when (plusp (nr-of-inner-dofs child-fe)) do
	  (loop for subcell across subcells and j from 0
		when (plusp (aref subcell-ndofs j)) do
		(loop with mblock = (make-real-matrix (nr-of-inner-dofs child-fe)
						      (aref subcell-ndofs j))
		      for child-dof in (fe-dofs child-fe)
		      when (interior-dof? child-dof) do
		      (loop for shape in (fe-basis fe) and dof in (fe-dofs fe)
			    when (= (dof-subcell-index dof) j) do
			    (setf (mref mblock (dof-in-vblock-index child-dof)
					   (dof-in-vblock-index dof))
				  (evaluate child-dof
					    (compose-2 shape (curry #'l2g child)))))
		      finally
		      (unless (mzerop mblock *ipt-tolerance*)
			(setf (tensor-ref imat i j) mblock)))))
    imat))

(defmethod compute-local-imatrix ((vecfe <vector-fe>))
  "Local imatrix computation for vector fe."
  (with-slots (refcell discretization components properties)
    vecfe
    (let* ((children (inner-refcell-children refcell))
	   (subcells (subcells refcell))
	   (subcell-ndofs (subcell-ndofs vecfe))
	   (subcell-offsets (subcell-offsets vecfe))
	   (children-offsets (children-offsets vecfe :inner))
	   (fe-imats (vector-map #'compute-local-imatrix components))
	   (vecfe-imat (make-instance '<sparse-tensor> :rank 2)))
      (loop for child across children and i from 0
	    for child-fe = (get-fe discretization child) do
	    (loop for subcell across subcells and j from 0
		  when (some #'(lambda (imat) (in-pattern-p imat i j))
			     fe-imats) do
		  (loop with mblock = (make-real-matrix (nr-of-inner-dofs child-fe)
							(aref subcell-ndofs j))
			for fe-imat across fe-imats and comp from 0
			for row-off = (aref (aref children-offsets comp) i)
			and col-off = (aref (aref subcell-offsets comp) j)
			when (in-pattern-p fe-imat i j) do
			(for-each-key-and-entry
			 #'(lambda (m n value)
			     (setf (mref mblock (+ m row-off) (+ n col-off))
				   value))
			 (tensor-ref fe-imat i j))
			finally (setf (tensor-ref vecfe-imat i j) mblock))))
      vecfe-imat)))

(defmemo local-imatrix (fe-disc refcell)
  "Memoized call of compute-local-imatrix."
  (compute-local-imatrix (get-fe fe-disc refcell)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Projection
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric compute-local-pmatrix (fe)
  (:documentation "Returns a local projection matrix as a sparse vector
with standard-matrix entries.  Note that projection is not as canonic as
interpolation.  We implement here the injection of dof functionals,
i.e. each degree of freedom is evaluated for the refined fe function."))

(defmethod compute-local-pmatrix ((fe <fe>))
  "This function computes a local projection matrix for Lagrangian finite
cells and maybe others.  For each interior degree of freedom, it finds the
child cell where the support of the node functional is contained and
evaluates it on the basis polynomials of the child which yields the entry
in the local projection matrix."
  (let ((nr-parent-inner-dofs (nr-of-inner-dofs fe)))
    (when (plusp nr-parent-inner-dofs)
      (let* ((refcell (reference-cell fe))
	     (fe-disc (discretization fe))
	     (all-children (refcell-children refcell))
	     (pmat (make-instance '<sparse-tensor> :rank 1)))
	(loop
	 for dof in (fe-dofs fe) and parent-dof-index below nr-parent-inner-dofs
	 for child = (find-if #'(lambda (child)
				  (and (= (dimension child) (dimension refcell))
				       (inside-cell? child (dof-coord dof))))
			      all-children) do
	 (loop
	  with child-fe = (get-fe fe-disc child)
	  and child-subcells = (subcells child)
	  for child-shape in (fe-basis child-fe)
	  and child-dof in (fe-dofs child-fe)
	  for entry = (evaluate dof (compose-2 child-shape (curry #'global->local child)))
	  do
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
		    entry)))))
	;; return pmat
	pmat))))

(defmethod compute-local-pmatrix ((vecfe <vector-fe>))
  "Computes a local projection matrix for vector-fe."
  (declare (optimize (debug 3)))
  (with-slots (refcell discretization components properties)
      vecfe
    (let ((nr-parent-inner-dofs (nr-of-inner-dofs vecfe)))
      (when (plusp nr-parent-inner-dofs)
	(let* ((all-children (refcell-children refcell))
	       (children-offsets (children-offsets vecfe :all))
	       (subcell-offsets (subcell-offsets vecfe))
	       (fe-pmats (vector-map #'compute-local-pmatrix components))
	       (vecfe-pmat (make-instance '<sparse-tensor> :rank 1)))
	  (loop for child across all-children and i from 0
	     for child-fe = (get-fe discretization child)
	     when (some #'(lambda (pmat) (and pmat (in-pattern-p pmat i))) fe-pmats) do
	     (loop with mblock = (make-real-matrix nr-parent-inner-dofs
						   (nr-of-inner-dofs child-fe))
		for pmat across fe-pmats and comp from 0
		for child-off = (aref (aref children-offsets comp) i)
		for parent-off = (aref (aref subcell-offsets comp) 0)
		when (and pmat (in-pattern-p pmat i)) do
		(for-each-key-and-entry
		 #'(lambda (m n value)
		     (setf (mref mblock (+ m parent-off) (+ n child-off))
			   value))
		 (tensor-ref pmat i))
		finally (setf (tensor-ref vecfe-pmat i) mblock)))
	  vecfe-pmat)))))

(defmemo local-pmatrix (fe-disc refcell)
  "Memoized call of compute-local-pmatrix."
  (compute-local-pmatrix (get-fe fe-disc refcell)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Transfer between different fe-spaces
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun local-transfer-matrix (fe-from fe-to)
  "Computes a local transfer matrix between different FE spaces."
  (assert (eq (reference-cell fe-from) (reference-cell fe-to)))
  (let* ((m (nr-of-dofs fe-to))
	 (n (nr-of-dofs fe-from))
	 (local-mat (make-real-matrix m n)))
    (loop for i from 0 and dof in (fe-dofs fe-to) do
	  (loop for j from 0 and phi in (fe-basis fe-from) do
		(setf (mref local-mat i j) (evaluate dof phi))))
    local-mat))
(memoize 'local-transfer-matrix)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing (test-fetransfer)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-fetransfer ()
  (let ((local-imat (compute-local-imatrix (get-fe (lagrange-fe 1) *unit-interval*))))
    (show local-imat)
    (tensor-ref local-imat 0 1))
  (show (compute-local-imatrix (get-fe (lagrange-fe 3) *unit-interval*)))
  (show (compute-local-pmatrix (get-fe (lagrange-fe 3) *unit-interval*)))
  (show (compute-local-imatrix (get-fe (lagrange-fe 2 :nr-comps 2) *unit-interval*)))
  (assert
   (= 28 (total-entries (compute-local-imatrix
			 (get-fe (lagrange-fe 2 :nr-comps 2) *unit-interval*)))))
  (show (compute-local-pmatrix (get-fe (lagrange-fe 4) *unit-interval*)))
  (show (compute-local-pmatrix (get-fe (lagrange-fe 1) *reference-vertex*)))
  (show (compute-local-pmatrix (get-fe (lagrange-fe 1 :nr-comps 2) *reference-vertex*)))
  (show (compute-local-pmatrix (get-fe (lagrange-fe 4 :nr-comps 2) *unit-interval*)))
  (show (compute-local-imatrix (get-fe (lagrange-fe 1) *unit-interval*)))
  (show (local-imatrix (lagrange-fe 1) *unit-interval*))
  (let ((count 0))
    (dotensor (entry (local-imatrix (lagrange-fe 1) *unit-interval*) :depth 2)
      (when entry (incf count)))
    (assert (= count 2)))
  )

;;; (test-fetransfer)
(fl.tests:adjoin-test 'test-fetransfer)
