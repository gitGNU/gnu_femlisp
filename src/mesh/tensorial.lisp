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

(defclass <tensorial> ()
  ()
  (:documentation "A mixin for simplex-product cells."))

(definline tensorial-p (obj) (typep obj '<tensorial>))

(defun tensorial-p (obj) (typep obj '<tensorial>))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Tensor product skeletons
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric make-product-cell (cell1 cell2 table)
  (:documentation "Generates the tensorial product of cell1 and cell2.  The
boundary is taken from the list of lower-dimensional products supplied in
table.  The tensor product boundary is in the order de1 x e2, e1 x de2."))

(defmethod make-product-cell :before ((cell1 <mapped-cell>) (cell2 <mapped-cell>) table)
  (error "Tensorial mappings do not yet work."))

(defgeneric tensor-product (x y)
  (:documentation
   "Generates the tensorial product of x and y.  x and y can be both
skeletons or cells."))

(defun product-table (skel1 skel2)
  (let ((product-table (make-hash-table :test #'equal)))
    (doskel (cell1 skel1 :direction :up)
      (doskel (cell2 skel2 :direction :up)
	(setf (gethash (list cell1 cell2) product-table)
	      (make-product-cell
	       cell1 cell2 product-table))))
    product-table))

(defmethod carthesian-product ((skel1 <skeleton>) (skel2 <skeleton>))
  ;; map product table to skeleton
  (loop with skel = (make-instance '<skeleton> :dimension
				   (+ (dimension skel1) (dimension skel2)))
	for cell being the hash-values of (product-table skel1 skel2)
	do (setf (skel-ref skel cell) nil)))

(defmethod carthesian-product ((cell1 <cell>) (cell2 <cell>))
  (gethash (list cell1 cell2)
	   (product-table (skeleton cell1) (skeleton cell2))))

;;; special cases: products with vertices
(defmethod make-product-cell ((vtx1 <vertex>) (vtx2 <vertex>) table)
  (declare (ignore table))
  (make-vertex (concatenate 'double-vec (vertex-position vtx1) (vertex-position vtx2))))

(defmethod make-product-cell ((vtx <vertex>) (cell <cell>) table)
  (make-instance
   (class-of cell)
   :boundary (map 'cell-vec #'(lambda (side) (gethash (list vtx side) table))
		  (boundary cell))))

(defmethod make-product-cell ((cell <cell>) (vtx <vertex>) table)
  (make-instance
   (class-of cell)
   :boundary (map 'cell-vec #'(lambda (side) (gethash (list side vtx) table))
		  (boundary cell))))

(defmethod make-product-cell ((cell1 <cell>) (cell2 <cell>) table)
  (make-instance
   (tensorial-class
    (mapcar #'dimension
	    (append (factor-simplices cell1) (factor-simplices cell2))))
   :boundary
   (concatenate 'cell-vec
		(map 'cell-vec #'(lambda (side) (gethash (list side cell2) table))
		     (boundary cell1))
		(map 'cell-vec #'(lambda (side) (gethash (list cell1 side) table))
		     (boundary cell2)))))

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

(defun weight-vector-factors (factors local-pos)
  "Given a local-pos this procedure computes a vector of weights which gives the
global position by taking the scalar product with the list of corners.  It
should be very useful to memoize the results of this function, because it is
probably called only on relatively few values of local-pos."
  (if (single? factors)
      (euclidean->barycentric local-pos)
      (let* ((dim1 (dimension (car factors)))
	     (weights-of-edges
	      (weight-vector-factors
	       (cdr factors)
	       (make-array (- (length local-pos) dim1) :element-type 'double-float
			   :displaced-to local-pos :displaced-index-offset dim1))))
	(apply #'concatenate 'double-vec
	       (map 'list #'(lambda (factor) (vec-s* factor weights-of-edges))
		    (euclidean->barycentric
		     (make-array dim1 :element-type 'double-float :displaced-to local-pos)))))))

(defun weight-vector-tensorial (cell local-pos)
  "Given a local-pos this procedure computes a vector of weights which gives the
global position by taking the scalar product with the list of corners.  It
should be very useful to memoize the results of this function, because it is
probably called only on relatively few values of local-pos."
  (declare (type <cell> cell) (type double-vec local-pos))
  (weight-vector-factors (factor-simplices cell) local-pos))

;;; the following speeds up bl-cdr by about 5%
;;; (memoize-symbol 'weight-vector-tensorial :test 'equalp)

(defmethod l2g ((cell <tensorial>) (local-pos array))
  (loop	with result = (make-double-vec (manifold-dimension cell))
	for weight across (weight-vector-tensorial (reference-cell cell) local-pos)
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
	collect (coerce (euclidean->barycentric 
			 (make-array (dimension factor) :element-type 'double-float
				     :displaced-to local-pos :displaced-index-offset pos))
			'list)))

(defun weight-lists-grad-simplex (simplex)
   "The result of this function are weight-lists with each weight-list
corresponding to the weights for a partial derivative."
  (mapcar
   #'(lambda (index)
       (mapcar
	#'(lambda (barycentric-index)
	    (cond ((zerop barycentric-index) -1.0d0)
		  ((= barycentric-index index) 1.0d0)
		  (t 0.0d0)))
	(range 0 (dimension simplex))))
   (range 1 (dimension simplex))))

(defun weight-lists-grad-tensorial (tensorial local-pos)
  (let* ((factors (factor-simplices tensorial))
	 (weight-lists (weight-lists factors local-pos)))
    (loop for factor in factors
	  and pos = 0 then (+ pos (dimension factor))
	  and weight-lists-tail on weight-lists
	  nconcing
	  (loop with temp = (car weight-lists-tail)
		for partial-derivative-wl in
		(weight-lists-grad-simplex factor)
		collecting
		(prog2 (shiftf temp (car weight-lists-tail) partial-derivative-wl)
		    (apply #'map-product #'* weight-lists)
		  (setf (car weight-lists-tail) temp))))))

(defmethod l2Dg ((cell <tensorial>) (local-pos array))
  (let* ((corners (corners cell))
	 (Dg (make-float-matrix (length (car corners)) (dimension cell))))
    (loop for weight-list in (weight-lists-grad-tensorial
			      (reference-cell cell) local-pos)
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

(defmethod generate-refine-info ((refcell <tensorial>))
  "Allocates refine-info vector with barycentric-corners, then fills the
boundary paths."
  (with-cell-information (refine-info)
    refcell
    (setq refine-info
	  (coerce
	   (apply #'map-product
		  #'(lambda (&rest child-factors)
		      (make-<child-info>
		       :reference-cell
		       (ensure-tensorial
			(remove 0 (mapcar (compose #'dimension #'child-reference-cell)
					  child-factors)))
		       :barycentric-corners
		       (apply #'append (mapcar #'child-barycentric-corners child-factors))
		       :boundary-paths ()))
		  (mapcar #'(lambda (factor) (coerce (refine-info factor) 'list))
			  (factor-simplices refcell)))
	   'child-info-vec)))
  (create-boundary-paths-of-tensorial refcell))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Tensorial class definition/generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-reference-tensorial (factor-dims)
  (reduce #'carthesian-product
	  (mapcar #'ensure-simplex factor-dims)))

(defun tensorial-class (factor-dims &optional mapped)
  (let* ((class-name (intern (format nil "<~{~S-~}TENSORIAL>" factor-dims) :mesh))
	 (class (find-class class-name nil)))
    (cond (class (if mapped (mapped-cell-class class) class))
	  (t
	   (let ((class (eval `(defclass ,class-name (<tensorial> <standard-cell>) ()))))
	      (let ((refcell (make-reference-tensorial factor-dims)))
		(initialize-cell-class
		 refcell (mapcar #'ensure-simplex factor-dims)))
	      class)))))

(defun ensure-tensorial (factor-dims)
  "Returns the reference tensorial for the given factor dimensions."
  (cond ((null factor-dims) *reference-vertex*)
	((single? factor-dims) (n-simplex (car factor-dims)))
	(t (reference-cell (tensorial-class factor-dims)))))

(defun n-cube (dim)
  "Returns the reference cube of dimension dim."
  (ensure-tensorial (make-list dim :initial-element 1)))

;;; immediate generation of commonly used tensorials, check for consistency
(defparameter *unit-quadrangle* (n-cube 2))
(defparameter *unit-cube* (n-cube 3))
(defparameter *unit-prism-1-2* (ensure-tensorial '(1 2)))
(defparameter *unit-prism-2-1* (ensure-tensorial '(2 1)))

;;;; Testing

(defun test-tensorial ()
  (n-cube 4)
  (ensure-tensorial '(3 1))
  (ensure-tensorial '(2 2))
  (ensure-tensorial '(1 3))
  (describe *unit-quadrangle*)
  (= 0.5 (aref (global->embedded-local (aref (boundary *unit-quadrangle*) 0) #(0.5 0.5)) 0))
  (describe (refcell-skeleton *unit-quadrangle*))
  (describe (refcell-refinement-skeleton *unit-quadrangle* 1))
  (refine-info *unit-cube*)
  (refcell-refinement-skeleton *unit-cube* 1)
  (describe (refcell-skeleton *unit-cube*))
  ;;
  (let ((child (aref (refine-info *unit-cube*) 8)))
    (format t "{D=~A} ~A :  ~A ~%"
	    (dimension (child-reference-cell child))
	    (child-barycentric-corners child)
	    (child-boundary-paths child)))
  (skeleton *unit-quadrangle*)
  (refine-globally (skeleton *unit-cube*))
  (global->embedded-local (aref (boundary *unit-quadrangle*) 0) #(0.5 0.5))
  (describe (refine-globally (skeleton *unit-cube*)))
  (refcell-refinement-vertices *unit-quadrangle* 2)
  )

;;; (test-tensorial)
(tests::adjoin-femlisp-test 'test-tensorial)