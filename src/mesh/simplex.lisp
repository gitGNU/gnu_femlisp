;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; simplex.lisp
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
;;;; class <simplex>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<simplex> (:conc-name simplex-) (:include <cell>)
		      (:print-function print-simplex))
  "The simplicial cell structure.")

#+ignore
(defclass <simplex> (<cell>)
  ()
  (:documentation "The simplicial cell structure."))

(defun print-simplex (simplex stream depth)
  (declare (ignore depth))
  (print-unreadable-object
   (simplex stream :type t :identity t)
   (when *print-cell*
     (whereas ((cell-class (cell-class simplex)))
       (format stream "{DIM=~A}" (cell-class-dimension cell-class)))
     (case *print-cell*
       (:corners (format stream "{MP=~A}" (midpoint simplex)))))))

(defun simplex? (obj)
  (or (typep obj '<simplex>)
      (typep obj '<vertex>)))

(defun reference-simplex (dim)
  (find-reference-cell-with
   #'(lambda (cell)
       (and (= (dimension cell) dim) (simplex? cell)))))

(defun simplex-class (dim)
  (whereas ((refcell (reference-simplex dim)))
    (cell-class refcell)))

(defun reference-simplices ()
  (find-reference-cells-with #'simplex?))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; vertices
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The following method requires a global consistent ordering of the simplex
;;; and its skeleton.  This makes the accesses to the vertices faster and seems
;;; to be necessary for the Freudenthal/Bey simplex refinement in more than
;;; three-dimensional problems).
(defmethod vertices ((simplex <simplex>))
  (let ((bdry (boundary simplex)))
    (cons (car (vertices (aref bdry 1)))
	  (vertices (aref bdry 0)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; local->global
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun euclidean->barycentric (pos)
  (let ((vec (make-double-vec (1+ (length pos)))))
    (setf (aref vec 0) (reduce #'- pos :initial-value 1.0d0))
    (dotimes (i (length pos) vec)
      (setf (aref vec (1+ i)) (aref pos i)))))

;;; the weight-vector of local-pos gives the weight for each corner for
;;; evaluation of the linear cell mapping at local-pos
(defmethod weight-vector ((simplex <simplex>) (local-pos array))
  (euclidean->barycentric local-pos))

;;; local->global for simplices

;;; This is probably a kludge:-(. We allow situation where boundaries may be
;;; parametrized but the cell itself is not.  In that case we take the global
;;; position of the corresponding vertices.

(defmethod l2g ((cell <simplex>) (local-pos array))
  (loop with result = (make-double-vec (manifold-dimension cell))
	for bc across (euclidean->barycentric local-pos)
	and corner in (corners cell) do
	(x+=s*y result bc corner)
	finally (return result)))

;;; only a first approximation...
(defmethod l2Dg ((simplex <simplex>) (local-pos array))
  (let* ((dim (dimension simplex))
	 (corners (corners simplex))
	 (origin (car corners))
	 (mdim (length origin))
	 (mat (make-float-matrix mdim dim)))
    (loop for corner in (cdr corners)
	  for col from 0 do
	  (dotimes (row mdim)
	    (setf (mat-ref mat row col)
		  (- (aref corner row) (aref origin row))))
	  finally (return mat))))

(defmethod l2jet ((simplex <simplex>) (local-pos array) (k integer))
  "The jet of a simplex with linear cell mapping is 0 above the first
derivative."
  (loop with dim = (dimension simplex)
	for i from 0 upto k
	collect
	(case i
	  ((0) (l2g simplex local-pos))
	  ((1) (l2Dg simplex local-pos))
	  (t (make-real-tensor (make-uint-vec i dim))))))

(defmethod coordinates-inside? ((cell <simplex>) local-pos)
  (and (not (some #'minusp local-pos))
       (<= (reduce #'+ local-pos) 1.0d0)))

(defmethod local-coordinates-of-midpoint ((cell <simplex>))
  (let ((dim (dimension cell)))
    (make-double-vec dim (/ 1.0d0 (1+ dim)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Regular simplex refinement by Freudenthal/Bey
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun two-sorted-parts (n)
  "Returns all partitions of {1,...,n} which consist of two sorted
parts.  This is used in Freudenthal/Bey's simplex refinement.
Example: (two-sorted-parts 2) -> ((() (1 2)) ((1) (2)) ((2) (1)) ((1 2) ()))"
  (loop with indices = (range 1 n)
	for k from 0 upto n
	nconc (loop for set in (k-subsets indices k)
		    collect (list set (ordered-set-difference indices set)))))

(defun freudenthal-refinement (n)
  "Returns all sub-simplices resulting from regular refinement of an
n-dimensional simplex by Freudenthal's algorithm.  The simplices are
given by their corners in barycentric coordinates (multiplied by 2 to
work with integer coordinates).
Example: (freudenthal-refinement 1) -> ((#(2 0) #(1 1)) (#(1 1) #(0 2)))"
  (let* ((n+1 (1+ n))
	 (positions (make-array (list n+1 n+1))))
    ;; Remark: we will use only one half (i<=j) of the positions-array to
    ;; allow maximal identification of objects.
    (for (i 0 n)
      (for (j i n)
	(setf (aref positions i j)
	      (vec-s* 0.5d0 (vec+ (unit-vector n+1 i) (unit-vector n+1 j))))))
    (loop for indices in (two-sorted-parts n)
	  collect
	  (let ((index-vector (permutation-shifted-inverse
			       (list->fixnum-vec (apply #'append indices))))
		(vi 0)
		(vj (length (car indices))))
	    (cons (aref positions vi vj)
		  (loop for l from 0 below n do
			(cond ((= vi (1- (aref index-vector l)))
			       (setq vi (aref index-vector l)))
			      ((= vj (1- (aref index-vector l)))
			       (setq vj (aref index-vector l)))
			      (t (error "something wrong in algorithm")))
			(if (> vi vj) (rotatef vi vj))
			collect (aref positions vi vj)))))))

(defun sub-cells-of-child (child-corners)
  (loop for i from 1 upto (length child-corners)
	collect (k-subsets child-corners i)))

(defun sub-cells-of-children (n)
  (loop with refinement = (make-list (1+ n))
	for child in (freudenthal-refinement n) do
	(setq refinement
	      (mapcar #'(lambda (set1 set2)
			  (ordered-union set1 set2 :test #'equalp))
		      refinement (sub-cells-of-child child)))
	finally (return refinement)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Generation of corresponding refinement rules
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
(defun get-path-create (simplex corners-of-sub-simplex path)
  "Returns the path (see the description in the definition of <child-info>) to
the sub-simplex given by its corners in barycentric coordinates.  Additionally,
if the sub-simplex cannot be found, it generates a new entry in the refine-info
vector and fills the barycentric-corners field."
  ;; if possible descend to a boundary
  (loop for side across (boundary simplex)
	and comp from 0
	when (every #'(lambda (vec) (zerop (aref vec comp)))
		    corners-of-sub-simplex)
	return (get-path-create
		side
		(mapcar #'(lambda (vec) (vector-cut vec comp))
			corners-of-sub-simplex)
		(cons comp path))
	finally  ; otherwise, search in the refine-info field
	(return
	  (loop for corners across (refine-info simplex)
		and i from 0
		when (equalp (car (child-barycentric-corners corners))
			     corners-of-sub-simplex)
		return (nreverse (cons i path))
		finally  ; the child was not found: create an appropriate entry
		(return
		  (let ((refine-info (refine-info simplex)))
		    (adjust-array refine-info (1+ (length refine-info))
				  :initial-element
				  (make-<child-info>
				   :class (cell-class
					   (reference-simplex
					    (1- (length corners-of-sub-simplex))))
				   :barycentric-corners (list corners-of-sub-simplex)
				   :boundary-paths nil))
		    ;; return reversed path
		    (nreverse (cons i path))))))))

(defun create-boundary-paths (simplex)
  "Create the boundary-paths for the <child-info> entry in the refine-info
vector.  We assume that the ordering of the simplex boundaries is opposite to
the ordering of the corners."
  (loop for child across (refine-info simplex) do
	(setf (child-boundary-paths child)
	      (let ((child-corners (car (child-barycentric-corners child))))
		(if (single? child-corners)
		    '()
		    (mapcar #'(lambda (corner)
				(get-path-create
				 simplex (remove corner child-corners :test #'equalp) '()))
			    child-corners))))))

(defmethod primary-refine-info ((refcell <simplex>))
  "We allocate an empty vector which we fill in update-refine-info."
  (make-array 0 :adjustable t))

(defmethod update-refine-info! ((refcell <simplex>))
  ;; fill refine-info with barycentric corners
  (loop for subcells-of-dim in (sub-cells-of-children (dimension refcell)) do
	(loop for corners-of-sub-simplex in subcells-of-dim do
	      (get-path-create refcell corners-of-sub-simplex '())))
  ;; ...and from there the boundary paths
  (create-boundary-paths refcell))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; simplex class definition/generation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun simplex-corner (dim i)
  "Creates the i-th corner of a simplex of dimension dim.  The 0-th corner is
the zero-vector.  The others are equal to (unit-vector dim 1-i)."
  (if (zerop i)
      (make-double-vec dim)
      (unit-vector dim (1- i))))

(defun make-reference-simplex (dim)
  "Constructs a reference simplex of the given dimension."
  (let ((corners (loop for i from 0 upto dim
		       collect (make-vertex (simplex-corner dim i))))
	(simplex-ht (make-hash-table :test #'equalp)))
    ;; intern the corners
    (loop for corner in corners do
	  (setf (gethash (list corner) simplex-ht) corner))
    ;; then intern the other cells
    (loop for subset in (k->l-subsets corners 2 (1+ dim)) do
	  (setf (gethash subset simplex-ht)
		(make-<simplex>
		 :cell-class (simplex-class (1- (length subset)))
		 :boundary (loop with bdry-vec = (make-array (length subset)
							     :element-type '<cell>
							     :initial-element *default-cell*)
				 for corner in subset and i from 0 do
				 (setf (aref bdry-vec i)
				       (gethash (remove corner subset :test #'equalp) simplex-ht))
				 finally (return bdry-vec)))))
    ;; the entry for corners is the desired simplex
    (gethash corners simplex-ht)))

(defun ensure-simplex (dim)
  "If the reference simplex of the given dimension exists, it is returned.
Otherwise, this function creates the simplex class, makes a reference simplex,
and initializes all its parameters.  All lower-dimensional simplices are
generated as well if it is necessary."
  (or (reference-simplex dim)
      (progn
	;; ensure that all lower dimensional simplices exist
	(for< (i 0 dim) (ensure-simplex i))
	;; generate reference cell and class
	(let* ((refcell (make-reference-simplex dim))
	       (cell-class (make-<cell-class> :factor-simplices (list refcell))))
	  (setf (cell-class refcell) cell-class)
	  ;; initialize
	  (initialize-cell-class refcell)
	  ;; finally, return the reference-cell
	  refcell))))

;;; access to commonly used simplices - this is also a check for consistency
(defparameter *unit-interval* (ensure-simplex 1))
(defparameter *1-simplex* (cell-class *unit-interval*))
(defparameter *unit-triangle* (ensure-simplex 2))
(defparameter *2-simplex* (cell-class *unit-triangle*))
(defparameter *unit-tetrahedron* (ensure-simplex 3))
(defparameter *3-simplex* (cell-class *unit-tetrahedron*))

;;; and lazy access to rarely used ones
(defun n-simplex (dim) (ensure-simplex dim))
(defun n-simplex-class (dim) (cell-class (ensure-simplex dim)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Mesh construction
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun make-simplex (boundary &key (check t) mapping)
  (let ((dim (1- (length boundary))))
    (when check
      (assert (every #'(lambda (side)
			 (= (dimension side) (1- dim)))
		     boundary))
      (assert (or (= dim 1) (closed? (skeleton boundary)))))
    (make-<simplex>
     :cell-class (simplex-class dim)
     :boundary boundary
     :mapping mapping)))

(defun make-line (from-vtx to-vtx &key (check t) mapping)
  (make-simplex (vector to-vtx from-vtx) :check check :mapping mapping))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun test-simplex ()
  (refine-info *unit-triangle*)
  (describe (refcell-skeleton *unit-triangle*))
  (describe (refcell-refinement-skeleton *unit-triangle* 1))
  (describe (refcell-refinement-skeleton *unit-triangle* 1))
  (refcell-refinement-skeleton *unit-triangle* 2)
  (mapcar #'corners
	  (cells-of-highest-dim
	   (refcell-refinement-skeleton *unit-triangle* 1)))
  (l2g *unit-triangle* #(0.5d0 0.5d0))
  (let ((bc (euclidean->barycentric #(0.5d0 0.5d0)))
	(corners (corners *unit-triangle*)))
    (mapc #'vec-s* (coerce bc 'list) corners))
  (local->global *unit-triangle* #(0.5d0 0.5d0))
  (n-simplex 5)
  (describe (refine-globally (skeleton *unit-tetrahedron*)))
  (refine-info *unit-interval*)
  )

;;; (test-simplex)
(tests::adjoin-femlisp-test 'test-simplex)

