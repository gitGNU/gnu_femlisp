;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  tensorial.lisp
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; class <tensorial>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<tensorial> (:conc-name tensorial-) (:include <cell>)
			(:print-function print-tensorial))
  "Class of a tensorial element.  Topologically, this is an arbitrary
product of simplices.")

#+ignore
(defclass <tensorial> (<cell>)
  ()
  (:documentation "Class of a tensorial element.  Topologically, this is an
arbitrary product of simplices."))

(defun print-tensorial (tensorial stream depth)
  (declare (ignore depth))
  (print-unreadable-object
   (tensorial stream :type t :identity t)
   (when *print-cell*
     (whereas ((cell-class (cell-class tensorial)))
       (format stream "{DIM=~A}"
	       (mapcar #'dimension (cell-class-factor-simplices cell-class))))
     (case *print-cell*
       (:corners (format stream "{MP=~A}" (midpoint tensorial)))))))

(defun tensorial? (obj) (typep obj '<tensorial>))
(defun cube? (obj)
  (and (typep obj '<cell>)
       (every #'(lambda (factor) (= (dimension factor) 1))
	      (factor-simplices obj))))

(defun reference-tensorial (dimensions)
  (find-reference-cell-with
   #'(lambda (cell)
       (equal (mapcar #'dimension (factor-simplices cell))
	      dimensions))))

(defun tensorial-class (dimensions)
  (aand (reference-tensorial dimensions)
	(cell-class it)))

(defun reference-tensorials ()
  (find-reference-cells-with #'tensorial?))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Tensor product skeletons
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric make-tensorial-cell (cell1 cell2 table)
  (:documentation "Generates the tensorial product of cell1 and cell2.  The
boundary is taken from the list of lower-dimensional products supplied in
table.  The tensor product boundary is in the order de1 x e2, e1 x de2."))

;;; an error check
(defmethod make-tensorial-cell
    :before ((cell1 <cell>) (cell2 <cell>) table)
    (declare (ignore table))
    (when (or (and (not (vertex? cell1)) (mapping cell1))
	      (and (not (vertex? cell2)) (mapping cell2)))
	(error "Tensorial mappings do not yet work.")))

(defgeneric tensor-product (x y)
  (:documentation
   "Generates the tensorial product of x and y.  x and y can be both
skeletons or cells."))

(defun tensor-product-table (skel1 skel2)
  (let ((product-table (make-hash-table :test #'equal)))
    (doskel (cell1 skel1 :direction :up)
      (doskel (cell2 skel2 :direction :up)
	(setf (gethash (list cell1 cell2) product-table)
	      (make-tensorial-cell
	       cell1 cell2 product-table))))
    product-table))

(defmethod tensor-product ((skel1 <skeleton>) (skel2 <skeleton>))
  ;; map product table to skeleton
  (loop with skel = (make-instance '<skeleton> :dimension
				   (+ (dimension skel1) (dimension skel2)))
	for cell being the hash-values of (tensor-product-table skel1 skel2)
	do (setf (skel-ref skel cell) nil)))

(defmethod tensor-product ((cell1 <cell>) (cell2 <cell>))
  (gethash (list cell1 cell2)
	   (tensor-product-table (skeleton cell1) (skeleton cell2))))

;;; special cases: products with vertices
(defmethod make-tensorial-cell ((vtx1 <vertex>) (vtx2 <vertex>) table)
  (declare (ignore table))
  (make-vertex (concatenate 'double-vec (vertex-position vtx1) (vertex-position vtx2))))

(defmethod make-tensorial-cell ((vtx <vertex>) (cell <cell>) table)
  (let ((new-cell (copy-structure cell)))
    (setf (boundary new-cell)
	  (map 'cell-vec #'(lambda (side) (gethash (list vtx side) table))
	       (boundary cell)))
    new-cell))

(defmethod make-tensorial-cell ((cell <cell>) (vtx <vertex>) table)
  (let ((new-cell (copy-structure cell)))
    (setf (boundary new-cell)
	  (map 'cell-vec #'(lambda (side) (gethash (list side vtx) table))
	       (boundary cell)))
    new-cell))

;;; general case
(defmethod make-tensorial-cell ((cell1 <cell>) (cell2 <cell>) table)
  (let ((factors (append (factor-simplices cell1) (factor-simplices cell2))))
    (make-<tensorial>
     :cell-class
     (whereas ((refcell (find-reference-cell-from-factors factors)))
       (cell-class refcell))
     :boundary
     (concatenate 'cell-vec
		  (map 'cell-vec #'(lambda (side) (gethash (list side cell2) table))
		       (boundary cell1))
		  (map 'cell-vec #'(lambda (side) (gethash (list cell1 side) table))
		       (boundary cell2))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vertices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Corners of tensorials are products of corners of their factors.  Therefore,
;;; a reasonable ordering of tensorial corners is the lexicographical ordering
;;; wrt the ordering of the corners of their factors.
(defmethod pseudo-vertices ((cell <cell>) (factor <vertex>))
  (list (vertices cell)))
(defmethod pseudo-vertices ((cell <tensorial>) (factor <simplex>))
  (let ((bdry (boundary cell))
	(factor-bdry (boundary factor)))
    (cons (car (pseudo-vertices (aref bdry 1) (aref factor-bdry 1)))
	  (pseudo-vertices (aref bdry 0) (aref factor-bdry 0)))))
(defmethod vertices ((cell <tensorial>))
  (apply #'append (pseudo-vertices cell (car (factor-simplices cell)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; l2g
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod weight-vector ((factors list) (local-pos array))
  "Given a local-pos this procedure computes a vector of weights which gives the
global position by taking the scalar product with the list of corners.  It
should be very useful to memoize the results of this function, because it is
probably called only on relatively few values of local-pos."
  (if (single? factors)
      (euclidean->barycentric local-pos)
      (let* ((dim1 (dimension (car factors)))
	     (weights-of-edges
	      (weight-vector (cdr factors)
			     (make-array (- (length local-pos) dim1) :element-type 'double-float
					 :displaced-to local-pos :displaced-index-offset dim1))))
	(apply #'concatenate 'double-vec
	       (map 'list #'(lambda (factor) (vec-s* factor weights-of-edges))
		    (euclidean->barycentric
		     (make-array dim1 :element-type 'double-float :displaced-to local-pos)))))))

(defmethod l2g ((cell <tensorial>) (local-pos array))
  (loop	with result = (make-double-vec (manifold-dimension cell))
	for weight across (weight-vector (factor-simplices cell) local-pos)
	and corner in (corners cell) do
	(x+=s*y result weight corner)
	finally (return result)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; l2Dg
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun weight-lists (simplices local-pos)
  "A list of weights for several simplices at once."
  (loop for factor in simplices
	and pos = 0 then (+ pos (dimension factor))
	collect (map 'list #'identity
		     (weight-vector
		      factor (make-array (dimension factor) :element-type 'double-float
					 :displaced-to local-pos :displaced-index-offset pos)))))

(defgeneric weight-lists-grad (cell local-pos)
  (:documentation "The result of this function are (dimension cell) weight-lists
with each weight-list corresponding to the weights for a partial
derivative."))

(defmethod weight-lists-grad ((simplex <simplex>) (local-pos array))
  (mapcar
   #'(lambda (index)
       (mapcar
	#'(lambda (barycentric-index)
	    (cond ((zerop barycentric-index) -1.0d0)
		  ((= barycentric-index index) 1.0d0)
		  (t 0.0d0)))
	(range 0 (dimension simplex))))
   (range 1 (dimension simplex))))

(defmethod weight-lists-grad ((tensorial <tensorial>) (local-pos array))
  (let* ((factors (factor-simplices tensorial))
	 (weight-lists (weight-lists factors local-pos)))
    (loop for factor in factors
	  and pos = 0 then (+ pos (dimension factor))
	  and weight-lists-tail on weight-lists
	  nconcing
	  (loop with temp = (car weight-lists-tail)
		for partial-derivative-wl in
		(weight-lists-grad
		 factor (make-array (dimension factor) :element-type 'double-float
				    :displaced-to local-pos :displaced-index-offset pos))
		collecting
		(prog2 (shiftf temp (car weight-lists-tail) partial-derivative-wl)
		    (apply #'map-product #'* weight-lists)
		  (setf (car weight-lists-tail) temp))))))

(defmethod l2Dg ((cell <tensorial>) (local-pos array))
  (let* ((corners (corners cell))
	 (Dg (make-float-matrix (length (car corners)) (dimension cell))))
    (loop for weight-list in (weight-lists-grad cell local-pos)
	  and i from 0 do
	  (loop for deriv across (reduce #'vec+ (mapcar #'vec-s* weight-list corners))
		and j from 0 do
		(setf (mat-ref Dg j i) deriv)))
    Dg))

(defmethod local-coordinates-of-midpoint ((cell <tensorial>))
  (apply #'concatenate 'double-vec
	 (mapcar #'local-coordinates-of-midpoint (factor-simplices cell))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; coordinates-inside?
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod coordinates-inside? ((cell <tensorial>) local-pos)
  (loop for factor in (factor-simplices cell)
	for factor-dim = (dimension factor)
	and offset = 0 then (+ offset factor-dim)
	unless (inside-cell? factor (vector-slice local-pos offset factor-dim))
	do (return nil)
	finally (return t)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Regular refinement of tensorials
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The generation of <refine-info> for tensorials assumes that the
;;; refinement data for reference simplices and lower-dimensional
;;; tensorials has already been generated.  The children of a tensorial are
;;; taken as tensorials of the factor children.  The barycentric-corners
;;; are a list of the barycentric corners in each factor.  Finally, the
;;; boundary-paths are computed.  This computation assumes that the
;;; boundaries of a tensorial appear in the order of their
;;; factor-simplices.

(defun create-boundary-paths-of-tensorial (cell)
  "Create the boundary-paths for the <child-info> entry in the
refine-info vector."
  (loop
   for child across (refine-info cell) do
   (setf (child-boundary-paths child)
	 (loop
	  for cell-factor in (factor-simplices cell)
	  for bc-factors on (child-barycentric-corners child)
	  for bc-factor = (car bc-factors)
	  unless (single? bc-factor)	; no boundary for a nodal factor
	  nconcing
	  (loop
	   for corner in bc-factor collect
	   (labels
	       ((find-path (cell cell-factor factor-side-corners)
		  (dbg :mesh "~&cell=~A~%  cell-factor=~A~%  factor-side-corners=~A~%"
		       cell cell-factor factor-side-corners)
		  (let ((path-to-bdry (and factor-side-corners
					   (get-path-create cell-factor factor-side-corners '()))))
		    (if (<= (length path-to-bdry) 1)
			;; in the interior
			(let ((pos (position
				    (append before-factors
					    (and factor-side-corners (list factor-side-corners))
					    (cdr bc-factors))
				    (refine-info cell) :key #'child-barycentric-corners :test #'equalp)))
			  (list pos))
			;; in the interior of a boundary face (bug: apparently not all covered!)
			(let* ((factor-side-id (car path-to-bdry))
			       (side-id (+ bdry-pos factor-side-id)))
			  (assert (< side-id (length (boundary cell))))
			  (cons
			   side-id
			   (let* ((factor-bcs-in-face
				   (unless (= (dimension cell-factor) 1)
				     (mapcar #'(lambda (vec) (vector-cut vec factor-side-id))
					     factor-side-corners)))
				  (factor-side-rc (reference-cell
						   (aref (boundary cell-factor) factor-side-id))))
			     (find-path (aref (boundary cell) side-id)
					factor-side-rc
					factor-bcs-in-face))))))))
	     (find-path cell cell-factor (remove corner bc-factor :test #'equalp))))
	  summing (nr-of-sides cell-factor) into bdry-pos
	  collecting bc-factor into before-factors))))

(defmethod primary-refine-info ((refcell <tensorial>))
  "Allocates refine-info vector with barycentric-corners."
  (list->vector
   'child-info-vec
   (apply #'map-product
	  #'(lambda (&rest child-factors)
	      (make-<child-info>
	       :class
	       (cell-class (find-reference-cell-from-factor-dimensions
			    (remove-if #'zerop
				       (mapcar #'(lambda (child)
						   (cell-class-dimension (child-class child)))
					       child-factors))))
	       :barycentric-corners
	       (apply #'append (mapcar #'child-barycentric-corners child-factors))
	       :boundary-paths ()))
	  (mapcar #'(lambda (factor) (vector->list (refine-info factor)))
		  (factor-simplices refcell)))))

(defmethod update-refine-info! ((refcell <tensorial>))
  "Fills boundary paths before doing the standard update for the
refine-info vector."
  (create-boundary-paths-of-tensorial refcell))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tensorial class definition/generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-reference-tensorial (factor-dims)
  (reduce #'make-tensorial-cell
	  (mapcar #'reference-simplex factor-dims)))

(defun make-tensorial-class-and-refcell (factor-dims)
  "This function creates the tensorial class, generates a reference cell,
and initializes all its parameters.  Lower-dimensional tensorials are also
generated if necessary."
  ;; ensure tensorials for boundary
  (loop for fdims on factor-dims
	for fdim = (car fdims) do
	(ensure-tensorial
	 (remove 0 (append before-fdims (list (1- fdim)) (cdr fdims))))
	collect fdim into before-fdims)
  ;; generate reference cell and class
  (let* ((refcell (reduce #'tensor-product
			  (mapcar #'reference-simplex factor-dims)))
	 (factors (mapcar #'reference-simplex factor-dims))
	 (cell-class (make-<cell-class> :factor-simplices factors)))
    (setf (cell-class refcell) cell-class)
    (initialize-cell-class refcell)
    refcell))

(defun ensure-tensorial (factor-dims)
  "If the reference tensorial exists, it is returned.  Otherwise, this function
calls make-tensorial-class-and-refcell to generate a new tensorial class and its
reference cell."
  (if (single? factor-dims)
      (ensure-simplex (car factor-dims))
      (or (find-reference-cell-from-factor-dimensions factor-dims)
	  (make-tensorial-class-and-refcell factor-dims))))

;;; immediate generation of commonly used tensorials, check for consistency
(defparameter *unit-quadrangle* (ensure-tensorial '(1 1)))
(defparameter *1-1-tensorial* (cell-class *unit-quadrangle*))
(defparameter *unit-cube* (ensure-tensorial '(1 1 1)))
(defparameter *1-1-1-tensorial* (cell-class *unit-cube*))
(defparameter *unit-prism-1-2* (ensure-tensorial '(1 2)))
(defparameter *unit-prism-2-1* (ensure-tensorial '(2 1)))

(defun n-cube (dim)
  "Returns the reference cube of dimension dim."
  (ensure-tensorial (make-list dim :initial-element 1)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Mesh construction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-tensorial (cell-class boundary &key (check t) mapping)
  "Tensorial constructor which does some checks."
  (assert (typep (cell-class-reference-cell cell-class) '<tensorial>))
  (let ((dim (cell-class-dimension cell-class)))
    (when check
      (assert (= (length boundary) (cell-class-nr-of-sides cell-class)))
      (assert (every #'(lambda (side)
			 (= (dimension side) (1- dim)))
		     boundary))
      (assert (or (= dim 1) (closed? (skeleton boundary)))))
    (make-<tensorial>
     :cell-class cell-class
     :boundary boundary
     :mapping mapping)))

(defun make-cube (boundary &key (check t) mapping)
  "Constructs a suitable n-cube with the given boundary."
  (let ((refcell (find-reference-cell-with
		  #'(lambda (refcell)
		      (and (or (eq refcell *unit-interval*)
			       (typep refcell '<tensorial>))
			   (= (nr-of-sides refcell) (length boundary)))))))
    (if (eq refcell *unit-interval*)
	(make-simplex boundary :check check :mapping mapping)
	(make-tensorial (cell-class refcell) boundary :check check :mapping mapping))))


;;;; Testing

(defun test-tensorial ()
  (n-cube 4)
  (ensure-tensorial '(3 1))
  (ensure-tensorial '(2 2))
  (ensure-tensorial '(1 3))
  (make-tensorial-class-and-refcell '(1 3))
  (= 0.5 (aref (global->embedded-local (aref (boundary *unit-quadrangle*) 0) #(0.5 0.5)) 0))
  (describe (refcell-skeleton *unit-quadrangle*))
  (describe (refcell-refinement-skeleton *unit-quadrangle* 1))
  (refine-info *unit-cube*)
  (refcell-refinement-skeleton *unit-cube* 1)
  (describe (refcell-skeleton *unit-cube*))
  ;;
  (let ((child (aref (refine-info *unit-cube*) 8)))
    (format t "{D=~A} ~A :  ~A ~%"
	    (mapcar #'dimension (cell-class-factor-simplices (child-class child)))
	    (child-barycentric-corners child)
	    (child-boundary-paths child)))
  (skeleton *unit-quadrangle*)
  (refine-globally (skeleton *unit-cube*))
  (global->embedded-local (aref (boundary *unit-quadrangle*) 0) #(0.5 0.5))
  (describe (refine-globally (skeleton *unit-cube*)))
  (refcell-refinement-vertices *unit-quadrangle* 2)
  (describe (make-cube (boundary *unit-interval*)))
  (describe (make-cube (boundary *unit-quadrangle*)))
  )

;;; (test-tensorial)
(tests::adjoin-femlisp-test 'test-tensorial)