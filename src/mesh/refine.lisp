;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; refine.lisp
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

(in-package :mesh)

;;;; This module provides the basic structure of refinement information and
;;;; for the regular refinement of skeletons.  The actual form of
;;;; refinement information depends on the cell type and is provided
;;;; separately in the modules (simplex) and (tensorial).

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; <child-info>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<child-info> (:conc-name child-))
  "This structure contains information about the children, their factor
simplices, and how to find their boundary sides.  The slots are:

class:  the child's cell-class

barycentric-corners:
   A list of the child's corners in barycentric coordinates

boundary-paths:
   For each side of the child this list contains a path in the form (i_1,
..., i_l, j): here i_1, ..., i_l are positions in subsequent boundary lists
starting from the boundary of the child's parent, and j is the position in
the refine-info vector of the boundary's parent.

transform-A, transform-b:  determine the transformation mapping for the child"
  
  (class nil :type (or <cell-class> null))
  (barycentric-corners () :type list)
  (boundary-paths () :type list)
  (transform-A nil)
  (transform-b nil))

(deftype child-info-vec () '(simple-array <child-info> (*)))

(defparameter *refinement-symbols* '(PARENT CHILDREN)
  "This list contains symbols used in the skeleton property list of
cells.")

(definline children (cell skeleton)
  (the (or null (simple-array <cell> (*)))
    (getf (skel-ref skeleton cell) 'CHILDREN)))

(definline (setf children) (child-vec cell skeleton)
  (setf (getf (skel-ref skeleton cell) 'CHILDREN)
	child-vec))

(definline parent (cell skeleton)
  (the (or null <cell>)
    (getf (skel-ref skeleton cell) 'PARENT)))

(definline (setf parent) (parent cell skeleton)
  (setf (getf (skel-ref skeleton cell) 'PARENT)
	parent))

(definline refined-p (cell skeleton)
  (children cell skeleton))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; refinement-skeleton generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; This is called inside the generation of refine-info!

(let ((refcell-refinement-table (make-hash-table :test 'equal)))
  (defun refcell-refinement-skeleton (refcell &optional (level 1))
    "This function is needed in its whole generality e.g. for plotting.
Calling it for large level arguments may cause a breakdown of the program,
thus the level argument is bounded conservatively."
    (assert (reference-cell? refcell))
    (assert (<= level 4))
    (or (gethash (cons refcell level) refcell-refinement-table)
	(let ((result
	       (cond
		 ((zerop level) (refcell-skeleton refcell))
		 (t (refine-globally (refcell-refinement-skeleton refcell (1- level)))))))
	  (setf (gethash (cons refcell level) refcell-refinement-table)
		result)))))

(defun subcell-children (cell skeleton)
  "This procedure gets the children of all subcells of cell."
  (loop with vec = (make-array (nr-of-subcell-children cell) :initial-element nil)
	and k = 0
	for subcell across (subcells cell) do
	(loop for child across (children subcell skeleton) do
	      (setf (aref vec k) child)
	      (incf k))
	finally (return vec)))

(definline inner-refcell-children (refcell)
  "Returns the children of refcell."
  (children refcell (refcell-skeleton refcell)))

(definline refcell-children (refcell)
  "Returns the children for refcell and subcells."
  (subcell-children refcell (refcell-skeleton refcell)))

(defmemo refcell-refinement-index-table (refcell level)
  "Returns a hash-table vertices->indices for refinement skeletons of
reference cells.  This is needed for plotting."
  (let ((index-table (make-hash-table))
	(index -1))
    (skel-for-each-cell #'(lambda (vtx) (setf (gethash vtx index-table) (incf index)))
			(refcell-refinement-skeleton refcell level)
			:dimension 0)
    index-table))

(defmemo refcell-refinement-vertices (refcell level)
  "Transforms refcell-refinement-index-table into a vector."
  (assert (reference-cell? refcell))
  (let* ((index-table (refcell-refinement-index-table refcell level))
	 (vertex-array (make-array (hash-table-count index-table))))
    (maphash #'(lambda (vtx index) (setf (aref vertex-array index) vtx))
	     index-table)
    vertex-array))

(definline refcell-nr-of-subcell-children (refcell)
  (nr-of-cells (refcell-refinement-skeleton refcell 1)))

(defun nr-of-subcell-children (cell)
  (refcell-nr-of-subcell-children (reference-cell cell)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; refine-info generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric primary-refine-info (refcell)
  (:documentation "Allocates a refine-info vector for refcell and does a
partial setup."))

(defgeneric update-refine-info! (refcell)
  (:documentation "Fills derived refine-info components."))

(let ((refine-info-table (make-hash-table :test 'eq)))
  (defun refcell-refine-info (refcell)
    "Because results are memoized, this function should only be called for
reference cells."
    (or (gethash refcell refine-info-table)
	(prog1 (setf (gethash refcell refine-info-table)
		     (primary-refine-info refcell))
	  (update-refine-info! refcell)))))

(defun refine-info (cell)
  "This function may be called for every cell."
  (declare (values (array <child-info> (*))))
  (refcell-refine-info (reference-cell cell)))

(defmethod primary-refine-info :before (refcell)
  (assert (reference-cell? refcell)))

(defmethod update-refine-info! :after (refcell)
  "Generates the transformation mappings."
  ;; we have to ensure that the refined skeleton exists.  This refinement
  ;; of the standard cells needs the incomplete refine-info.
  (refcell-refinement-skeleton refcell 1)
  (loop for child-info across (refcell-refine-info refcell)
	and child across (inner-refcell-children refcell) do
	(let ((corner (car (corners child))))
	  (setf (child-transform-b child-info) corner)
	  (setf (child-transform-A child-info)
		(and (plusp (dimension child))
		     (l2Dg child (make-double-vec (dimension child))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Refinement of skeletons
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric refine-cell! (cell skel refined-skel refined-region)
  (:documentation "This local refinement function is used to assemble a
global refinement mapping which is represented by the skeleton
`refinement'.  It needs and ensures that the boundary of `cell' is
already refined.  An existing refinement of `cell' is simply kept."))

(defmethod refine-cell! ((vtx <vertex>) (skel <skeleton>) (refined-skel <skeleton>) refined-region)
  "This method may be modified by more specific methods.  Especially for
meshes the domain patches have to be put in the properties."
  (declare (optimize (speed 3)))
  (unless (refined-p vtx skel)
    (let* ((child (make-vertex (vertex-position vtx)))
	   (children-vector (vector child)))
      (when refined-region
	(setf (skel-ref refined-region vtx) children-vector))
      (setf (parent child refined-skel) vtx)
      (setf (children vtx skel) children-vector)
      nil)))

(defmethod refine-cell! ((cell <cell>) (skel <skeleton>) (refined-skel <skeleton>) refined-region)
  ;;(declare (optimize (speed 3)))
  (unless (refined-p cell skel)
    ;; first ensure that the boundary is already refined
    (loop for side across (boundary cell)
	  unless (refined-p side skel) do
	  (refine-cell! side skel refined-skel refined-region))
    ;; then allocate the children vector
    (let* ((refine-info (refine-info cell))
	   (nr-of-children (length refine-info))
	   (children-vector (make-array nr-of-children :initial-element nil)))
      ;; put the pair cell/children-vector already in the refined-skel
      (when refined-region
	(setf (skel-ref refined-region cell) children-vector))
      (setf (children cell skel) children-vector)
      ;; and fill this vector in place
      (dotimes (i nr-of-children children-vector)
	(setf (aref children-vector i)
	      (let* ((child-info (aref refine-info i))
		     (child-class (child-class child-info)))
		(if (eq child-class *vertex-class*)
		    ;; inner vertices appear only once in refinements of products of 1-simplices
		    (make-vertex (local->global cell (make-double-vec (dimension cell) 0.5d0)))
		    ;; we want to keep the class of the cell,
		    ;; therefore we copy a reference cell and
		    ;; reinitialize the boundary and mapping slot
		    (let ((new-cell (copy-structure
				     (cell-class-reference-cell child-class))))
		      (setf (boundary new-cell) 
			    (labels ((find-side (cell path)
				       (if (single? path)
					   (aref (children cell skel) (car path))
					   (find-side (aref (boundary cell) (car path))
						      (cdr path)))))
			      (map 'cell-vec #'(lambda (path) (find-side cell path))
				   (child-boundary-paths child-info))))
		      (setf (mapping new-cell)
			    (whereas ((cell-map (mapping cell)))
			      (transform-function
			       cell-map
			       (child-transform-A child-info)
			       (child-transform-b child-info))))
		      new-cell)))))
      ;; finally, insert the children in the refined skeleton.
      (loop for child across children-vector do
	    (setf (parent child refined-skel) cell)))))

(defmethod refine-cell! :after ((cell <cell>) (skel <skeleton>) (refined-skel <skeleton>) refined-region)
  "This after method handles the case where all identified cells have to be
refined.  All identified cells are refined if necessary and a new
identification is generated."
  (let ((identified-cells (cell-identification cell skel)))
    (when identified-cells
      ;; ensure refinement of identified cells; this is a recursive call
      (when (loop with every-p = t
		  for id-cell in identified-cells
		  unless (refined-p id-cell skel) do
		  (setq every-p nil)
		  (refine-cell! id-cell skel refined-skel refined-region)
		  finally (return every-p))
	;; all identified cells are refined, now set identification for
	;; all children
	(loop for i from 0 below (length (children cell skel))
	      for identified-children =
	      (loop for id-cell in identified-cells
		    collect (aref (children id-cell skel) i))
	      do
	      (loop for id-cell in identified-cells
		    for child = (aref (children id-cell skel) i) do
		    (setf (cell-identification child refined-skel) identified-children)))))))
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; global refinement
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; A refinement is a mapping from a skeleton to arrays of children.
;;; It is based on the cell refinement which assumes that the
;;; cell boundary is already refined.

(defun refine-globally (skel)
  "The refinement algorithm works as follows: It proceeds from
zero-dimensional vertices to higher-dimensional cells.  On each level of
the skeleton (starting from 0) the method refine-cell! is called on each
cell thus filling a refined-skel."
  (let ((refined-skel (make-analog skel)))
    (doskel (cell skel :direction :up)
      (refine-cell! cell skel refined-skel nil))
    refined-skel))

;;; For simplices and tensorial cells the construction of the refine-info
;;; and the refinement entry is done in the files simplex.lisp and
;;; tensorial.lisp.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; cell-class activation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; we add refcell-refine-info to the cell-class activation procedure
(adjoin-cell-class-initialization-function 'refcell-refine-info)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-skeleton-refinement ()
  "More tests are done when initializing other classes than vertex."
  (describe (refine-globally (skeleton *reference-vertex*)))
  )

;;; (test-skeleton-refinement)
(tests::adjoin-femlisp-test 'test-skeleton-refinement)
