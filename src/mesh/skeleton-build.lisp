;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; skeleton-build.lisp
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

(in-package :fl.mesh)

(defun skeleton-without-cell (skel cell-to-remove)
  "Removes a cell from a skeleton such that the rest remains a skeleton.
Warning: does not handle identifications yet."
  (let ((subcells-to-remove (subcells cell-to-remove)))
    (skeleton
     (find-cells #'(lambda (cell)
		     (let ((subcells (subcells cell)))
		       (not (or (find cell subcells-to-remove)
				(find cell-to-remove subcells)))))
		 skel))))

(defmethod copy-skeleton (skel &key properties transformation)
  "Copies a skeleton.  Properties is a list of properties to be copied."
  (let ((new-skel (make-analog skel))
	(table (make-hash-table)))
    (doskel (cell skel :direction :up)
      (let ((new-cell
	     (if (vertex? cell)
		 (make-vertex (copy-seq (vertex-position cell)))
		 (let ((copy (copy-cell cell)))
		   (setf (slot-value copy 'boundary)
			 (map 'cell-vec (rcurry #'gethash table) (boundary cell)))
		   copy))))
	(setf (skel-ref new-skel new-cell) ())
	(dolist (prop properties)
	  (whereas ((prop-val (getf (skel-ref skel cell) prop)))
	    (setf (getf (skel-ref new-skel new-cell) prop)
		  prop-val)))
	(setf (gethash cell table) new-cell)))
    ;; transfer identification information
    (doskel (cell1 skel)
    (let ((identified-cells1 (cell-identification cell1 skel)))
      (when identified-cells1
	(let ((cell2 (gethash cell1 table)))
	    (unless (cell-identification cell2 new-skel)
	      (let ((identified-cells2 (mapcar (rcurry #'gethash table)
					       identified-cells1)))
		(loop for cell2 in identified-cells2 do
		      (setf (cell-identification cell2 new-skel)
			    identified-cells2))))))))
    ;; when a transformation is defined, call it
    (when transformation
      (maphash transformation table))
    ;; return copy (and copy-table for further use)
    (values new-skel table)))

(defun refined-skeleton-copy (skel &optional (refinements 0))
  (loop for newskel = (copy-skeleton skel) then (refine-globally newskel)
	repeat refinements
	finally (return newskel)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Mesh construction in the UG way
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun insert-cell-from-corners (mesh corners->cell cell-class corners properties)
  "Creates a cell of type cell-class with corners given by corners.
corners->cell has to be an equalp hash-table mapping corners to the
corresponding cell.  It is updated by this function."
  (or (gethash corners corners->cell)
      (let* ((refcell (reference-cell cell-class))
	     (refcell-corners (corners refcell)))
	(assert (= (length refcell-corners) (length corners)))
	(flet ((refcell-corner->corner (refcell-corner)
		 (nth (position refcell-corner refcell-corners :test #'equal) corners)))
	  (let ((cell (if (vertex? refcell)
			  (make-vertex (car corners))
			  (copy-cell refcell))))
	    (unless (vertex-p cell)
	      (setf (slot-value cell 'boundary)
		    (vector-map
		     #'(lambda (refcell-side)
			 (let ((side-corners (mapcar #'refcell-corner->corner
						     (corners refcell-side))))
			   (insert-cell-from-corners
			    mesh corners->cell
			    (class-of refcell-side) side-corners properties)))
		     (boundary refcell))))
	    (setf (skel-ref mesh cell) properties)
	    (setf (gethash corners corners->cell) cell))))))

(defun structured-skeleton (N h &key corners->cell)
  "Create a uniform box skeleton consisting of N_1 x ... x N_dim cubes of
dimensions h_1 x ... x h_dim."
  (assert (= (length N) (length h)))
  (let* ((dim (length N))
	 (skel (make-instance '<skeleton> :dimension dim))
	 (cube-class (class-of (n-cube dim))))
    (ensure corners->cell (make-hash-table :test 'equalp))
    (multi-for (ivec (make-fixnum-vec dim 1) N)
      (let ((corners ()))
	(multi-for (jvec (make-fixnum-vec dim -1) (make-fixnum-vec dim 0))
	  (push (map 'double-vec
		     #'(lambda (k_i h_i) (float (* k_i h_i) 1.0))
		     (m+ ivec jvec) h)
		corners))
	(insert-cell-from-corners skel corners->cell cube-class (nreverse corners) ())))
    skel))

(defun skel-add! (skel-1 skel-2 &key (override ()) active-skel-1)
  "Adds skel-2 to skel-1 destructively for skel-1.  Overlaying objects are
identified.  Override is a list of properties which are taken from skel-2
on the overlap.  Active-skel-1 is used for hierarchical-meshes for
selecting a level to which skel-2 is added.  This function returns three
values: the first is skel-1, the second is the copy made of skel-2, the
third is a hash-table mapping overlapping cells from the copy of skel-2 to
their counterpart in skel-1."
  (let ((skel-2 (copy-skeleton skel-2 :properties override))
	(overlap (make-hash-table)))
    (doskel (cell-2 skel-2 :direction :up)
      (let* ((bdry (vector-map #'(lambda (side) (or (gethash side overlap) side))
			       (boundary cell-2)))
	     (twins
	      (find-cells
	       #'(lambda (cell-1)
		   (if (vertex? cell-1)
		       (equalp (vertex-position cell-1) (vertex-position cell-2))
		       (equalp (boundary cell-1) bdry)))
	       (or active-skel-1 skel-1) :dimension (dimension cell-2))))
	(cond ((null twins)		; add cell-2 to skel-1
	       (unless (vertex-p cell-2)
		 (setf (slot-value cell-2 'boundary) bdry))
	       (setf (skel-ref skel-1 cell-2) (skel-ref skel-2 cell-2)))
	      (t			; insert cell-2 in overlap
	       (unless (= (length twins) 1)
		 (error "skel-1 was already not overlap-free."))
	       (dolist (prop override)
		 (setf (get-cell-property (car twins) skel-1 prop)
		       (get-cell-property cell-2 skel-2 prop)))
	       (setf (gethash cell-2 overlap) (car twins))))))
    (values skel-1 skel-2 overlap)))

(defmethod transformed-skeleton ((skel <skeleton>) &key transformation properties)
  "Transforms skel by transforming the vertex positions."
  (copy-skeleton
   skel :properties properties :transformation
   #'(lambda (old-cell new-cell)
       (if (zerop (dimension old-cell))
	   (copy! (evaluate transformation (vertex-position old-cell))
		  (vertex-position new-cell))
	   (let ((mapping (cell-mapping old-cell)))
	     (change-class new-cell (mapped-cell-class (class-of new-cell))
			   :mapping (compose-2 transformation mapping)))))))

(defmethod linearly-transformed-skeleton ((skel <skeleton>) &key A b properties)
  "Transforms skel by transforming the vertex positions."
  (copy-skeleton
   skel :properties properties :transformation
   #'(lambda (old-cell new-cell)
       (if (zerop (dimension old-cell))
	   (copy! (m+ (m* A (vertex-position old-cell)) b)
		  (vertex-position new-cell))
	   (when (mapped-p old-cell)
	     (setf (mapping new-cell)
		   (transform-function
		    (mapping old-cell) :image-transform (list A b))))))))

(defmethod shift-skeleton ((skel <skeleton>) shift &key properties)
  "Shifts skel by vec.  vec has to be a vector of dimension
\(manifold-dimension skel\)."
  (linearly-transformed-skeleton
   skel :A (eye (dimension skel)) :b shift :properties properties))

(defmethod subskeleton ((skel <skeleton>) test)
  (let ((subskel (make-instance '<skeleton> :cells (find-cells test skel))))
    (doskel (cell subskel)
      (setf (skel-ref subskel cell) (skel-ref skel cell)))
    subskel))

(defun telescope (left-skel left->right)
  ;; generate products of left-skel cells with lines and set their boundaries
  (let ((product-table (make-hash-table :test #'equal))
	(left-node (aref (boundary *unit-interval*) 1))
	(right-node (aref (boundary *unit-interval*) 0)))
    ;; set left and right boundary
    (doskel (cell left-skel :direction :up)
      (setf (gethash (list left-node cell) product-table) cell)
      (setf (gethash (list right-node cell) product-table)
	    (gethash cell left->right)))
    ;; fill the space between the two skeleta
    (doskel (cell left-skel :direction :up)
      (let ((new-cell (make-product-cell *unit-interval* cell product-table)))
	(change-class new-cell
		      (mapped-cell-class (class-of new-cell))
		      :mapping
		      (homotopy (cell-mapping cell)
				(cell-mapping (gethash cell left->right))))
	(setf (gethash (list *unit-interval* cell) product-table)
	      new-cell)))
    ;; and generate the telescope skeleton
    (skeleton (hash-table-values product-table))))

;;; Testing: (test-skeleton-build)
(defun test-skeleton-build ()
  (describe
   (let ((skel (skeleton  *unit-quadrangle*)))
     (skel-add! (shift-skeleton skel #d(0.0 1.0)) skel)))
  (subskeleton (skeleton *unit-quadrangle*)
	       #'(lambda (cell) (= (aref (midpoint cell) 0) 0.0)))
)

(fl.tests:adjoin-test 'test-skeleton-build)
