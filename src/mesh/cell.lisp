;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;  cell.lisp
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; The standard cell subclass construction works as follows:
;;;;
;;;; 1. Generate the subclass.
;;;; 2. Generate a reference cell (with factor-simplices).
;;;; 3. Initialize the reference cell (puts it also in dictionary).
;;;; 4. The refinement information is generated as it is needed.
;;;;
;;;; This is done separately in vertex.lisp, simplex.lisp and
;;;; tensorial.lisp.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <cell> class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<cell-class> (:conc-name cell-class-))
  (dimension -1 :type fixnum)
  (nr-of-sides -1 :type fixnum)
  (nr-of-vertices -1 :type fixnum)
  ;; ordered subcell access
  (nr-of-subcells -1 :type fixnum)
  (subcell-parent-indices (fixnum-vec) :type fixnum-vec)
  (subcell-boundary-indices (fixnum-vec) :type fixnum-vec)
  (factor-simplices ())
  (reference-cell () :type t))

(deftype cell-vec () '(simple-array <cell> (*)))

(defstruct (<cell> (:conc-name nil) (:print-function print-<cell>))
  (cell-class nil :type (or <cell-class> null))
  (boundary #() :type cell-vec)
  (mapping nil))

#+ignore
(defclass <cell> ()
  ((cell-class :reader cell-class :initarg :cell-class :type (or <cell-class> null))
   (boundary :reader boundary :initarg :boundary :initform #() :type cell-vec)
   (mapping :reader mapping :initarg :mapping nil))
  (:documentation "The basic cell class.  Every cell consists of its class
collecting all class information, an array of boundary cells and a mapping
slot.  That slot contains the position for vertices and a possibly
nonlinear mapping for other cells.  A value of nil means that multilinear
interpolation between the corners is used for constructing the mapping."))

(defvar *default-cell* (make-<cell>))

(defparameter *print-cell* :corners
  "If set to :corners, prints the corners of a cell.")

(defun print-<cell> (cell stream depth)
  (declare (ignore depth))
  (print-unreadable-object
   (cell stream :type t :identity t)))

(defun copy-cell (cell)
  "Copy constructor for cells."
  (copy-structure cell))

(defmethod check ((cell <cell>))
  (loop with side-dim = (1- (dimension cell))
	for side across (boundary cell)	do
	(assert (= side-dim (dimension side))))
  (whereas ((mapping (and (not (vertex? cell)) (mapping cell))))
    ;; check if mapping is reasonable
    (when (differentiable-p mapping)
      (let* ((mp (local-coordinates-of-midpoint cell))
	     (grad (evaluate-gradient mapping mp))
	     (numgrad (evaluate (numerical-gradient mapping) mp)))
	(assert (< (norm (m- grad numgrad))
		   (* 1.0e-2 (max (norm grad) (norm numgrad)))))))))

;;; accessors
(declaim (ftype (function (*) positive-fixnum) dimension))
(defmethod dimension ((cell <cell>))
  (cell-class-dimension (cell-class cell)))
(defmethod nr-of-sides ((cell <cell>))
  (cell-class-nr-of-sides (cell-class cell)))
(defmethod nr-of-vertices ((cell <cell>))
  (cell-class-nr-of-vertices (cell-class cell)))
(defmethod nr-of-subcells ((cell <cell>))
  (cell-class-nr-of-subcells (cell-class cell)))
(defmethod factor-simplices ((cell <cell>))
  (cell-class-factor-simplices (cell-class cell)))
(defmethod factor-dimensions ((cell <cell>))
  (mapcar #'dimension (factor-simplices cell)))
(defmethod reference-cell ((cell <cell>))
  (cell-class-reference-cell (cell-class cell)))

(defmethod describe-object ((cell <cell>) stream)
  (format stream "~&~A~%" cell))

(defmethod describe-all ((cell <cell>))
  (describe cell)
  (mapc #'describe-all (boundary cell)))

; (defun print-cell (cell stream depth)
;   (print-unreadable-object (cell-class cell))
;   (print (boundary cell))
;   (print-unreadable-object (mapping cell)))

;;; Ordered access (i.e. the ordering depends only on the cells class) to
;;; sub-cells

(defun subcells (cell)
  "Returns a vector containing all subcells of a given cell."
  (declare (type <cell> cell))
  (let ((subcells (make-array (nr-of-subcells cell)))
	(parent-indices (cell-class-subcell-parent-indices (cell-class cell)))
	(boundary-indices (cell-class-subcell-boundary-indices (cell-class cell))))
    (setf (aref subcells 0) cell)
    (do ((k 1 (1+ k)))
	((= k (nr-of-subcells cell)) subcells)
      (setf (aref subcells k)
	    (aref (boundary (aref subcells (aref parent-indices k)))
		  (aref boundary-indices k))))))


(defun create-subcell-info (cell)
  "Initializes the slots necessary for subcell access."
  (let* ((all-subcells (list cell))
	 (subcell-tail all-subcells)
	 (parent-indices (list -1))
	 (boundary-indices (list -1))
	 (nr-of-subcells 1))
    (do ((subcells all-subcells (cdr subcells))
	 (parent-index 0 (1+ parent-index)))
	
	 ;; end of loops: set subcell info in class slots
	((null subcells)
	 (let ((cell-class (cell-class cell)))
	   (setf (cell-class-nr-of-subcells cell-class) nr-of-subcells)
	   (setf (cell-class-subcell-parent-indices cell-class)
		 (map 'fixnum-vec #'identity
		      (nreverse parent-indices)))
	   (setf (cell-class-subcell-boundary-indices cell-class)
		 (map 'fixnum-vec #'identity
		      (nreverse boundary-indices)))
	   cell))
      
      ;; boundary loop
      (dotimes (boundary-index (length (boundary (car subcells))))
	(let ((side (aref (boundary (car subcells)) boundary-index)))
	  (unless (member side all-subcells)
	    (incf nr-of-subcells)
	    (push parent-index parent-indices)
	    (push boundary-index boundary-indices)
	    (nconc subcell-tail (list side))))))))

(defmethod subcells-of ((cell <cell>) (dim fixnum))
  (labels ((test (subcell) (= dim (dimension subcell))))
    (let* ((subcells (subcells cell))
	   (first-with-dim (position-if #'test subcells))
	   (last-with-dim (position-if #'test subcells :from-end t)))
    (make-array (1+ (- last-with-dim first-with-dim))
		:displaced-to subcells
		:displaced-index-offset first-with-dim))))

;;; This could be defined in this way, but we want a special ordering in
;;; simplex.lisp and tensorial.lisp.  Thus, this method is overloaded in
;;; those files.
(defmethod vertices ((cell <cell>))
  (subcells-of cell 0)
  (error "This should be overloaded (because the cell mappings need another
ordering)."))

(defmethod corners ((cell <cell>))
  (mapcar #'(lambda (vtx) (local->global vtx (double-vec)))
	  (vertices cell)))

(defmethod g-corners ((cell <cell>))
  (mapcar #'(lambda (vtx) (l2g vtx (double-vec)))
	  (vertices cell)))

(defun is-cell? (cell corners &key (tolerance 0.0))
  "This routine is quite useful for debugging."
  (and (= (dimension cell) (1-(length corners)))
       (set-equal? corners (corners cell)
		   :test #'(lambda (xvec yvec)
			     (every #'(lambda (x y) (<= (abs (- x y)) tolerance))
				    xvec yvec)))))

(defmethod manifold-dimension ((cell <cell>))
  "Manifold-dimension is recursively defined, anchored at the definition
for vertices."
  (manifold-dimension (aref (boundary cell) 0)))

(defun cell-mapping (cell)
  "Returns either the mapping slot or, if it is NIL, a <special-function>
which is equivalent to calling l2g and l2Dg."
  (let ((mapping (mapping cell)))
    (if mapping
	(if (vectorp mapping)
	    (make-instance '<constant-function> :domain-dimension 0
			   :image-dimension (length mapping) :value mapping)
	    mapping)
	(make-instance
	 '<special-function>
	 :domain-dimension (dimension cell) :image-dimension (manifold-dimension cell)
	 :evaluator (curry #'l2g cell) :gradient (curry #'l2Dg cell)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; local->global for isoparametric cells (as an after method)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric local->global (cell local-pos)
  (:documentation "local->global checks if a mapping is given for the cell.
If yes, then this mapping is evaluated.  If no, then the function l2g is called
which should do a multilinear interpolation from the cell's corners."))

(defgeneric local->Dglobal (cell local-pos)
  (:documentation "local->Dglobal checks if a mapping is given for the cell.
If yes, then the gradient of this mapping is evaluated (if available).  If no,
then the function l2Dg is called which gives the gradient for a multilinear
interpolation from the cell's corners."))

(defmethod local->global ((cell <cell>) local-pos)
  (aif (mapping cell)
       (evaluate it local-pos)
       (l2g cell local-pos)))

(defmethod local->Dglobal ((cell <cell>) local-pos)
  (aif (mapping cell)
       (evaluate-gradient it local-pos)
       (l2Dg cell local-pos)))

;;; Only these functions (specialized to multilinear mappings) are
;;; different for each cell class.

(defgeneric l2g (cell local-pos)
  (:documentation "Computes the global position by interpolation from
the vertices."))
  
(defgeneric l2Dg (cell local-pos)
  (:documentation "Computes the gradient for a multilinear
interpolation from the vertices."))

;;; It may be faster to evaluate function/derivatives in one sweep.  The
;;; following interface is a trial version to provide this functionality.  For
;;; higher order problems and Hermite ansatz functions, it will probably be
;;; necessary to generalize this to return the k-jet of a function.

(defgeneric l2jet (cell coords k)
  (:documentation "Returns the k-jet of the multilinear interpolation from
vertices, i.e. multiple values g, Dg, ..., D^k g."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; global->local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric global->local (cell global-pos)
  (:documentation "Mainly useful for fe evaluation: from the local position, the
value of a fe function can be obtained by interpolation.  This is done by a
Newton iteration, which converges in one step for linear mappings."))

(defmethod global->local ((cell <cell>) global-pos)
  "Does a Newton iteration or a Gauss-Newton method for approximating
global-pos by the cell mapping."
  (if (= (dimension cell) (length global-pos))
      (newton #'(lambda (x) (m- (local->global cell x) global-pos))
	      (make-double-vec (dimension cell))
	      :approximate-gradient #'(lambda (x) (local->Dglobal cell x)))
      (let* ((A (local->Dglobal cell (make-double-vec (dimension cell))))
	     (At (transpose A))
	     (At*A (m* At A)))
	(newton #'(lambda (x) (m* At (m- (local->global cell x) global-pos)))
		(make-double-vec (dimension cell))
		:approximate-gradient (constantly At*A)))))

(defgeneric inside-cell? (cell global-pos)
  (:documentation "Checks if global-pos is inside the interior of the cell.
It calls coordinates-inside? which is defined for every cell class."))

(defgeneric coordinates-inside? (cell local-pos)
  (:documentation "Checks if the given local coordinates are inside the
reference cell."))

(defmethod inside-cell? ((cell <cell>) position)
  (let ((local (global->local cell position)))
    (coordinates-inside? cell local)))

(defgeneric local-coordinates-of-midpoint (cell)
  (:documentation "Returns local coordinates of the cell midpoint."))

(defun midpoint (cell)
  "Returns cell midpoint in global coordinates."
  (local->global cell (local-coordinates-of-midpoint cell)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; global->embedded-local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric global->embedded-local (cell global-pos)
  (:documentation "This computes a local coordinate which solves the
Ausgleichsproblem of mapping to a point as near as possible to global-pos.  It
involves more computational work than global->local.  As a second value, the
distance to global-pos is returned."))

(defmethod global->embedded-local ((cell <cell>) global-pos)
  (let ((x (newton #'(lambda (x)
		       (m* (transpose (local->Dglobal cell x))
			   (m- (local->global cell x) global-pos)))
		   (make-double-vec (dimension cell))
		   :approximate-gradient
		   #'(lambda (x)
		       (let ((A (local->Dglobal cell x)))
			 (m* (transpose A) A)))
		   :maxsteps 1)))
    (values x (norm (m- (local->global cell x) global-pos)))))

(defgeneric inside-cell? (cell global-pos)
  (:documentation "Checks if global-pos is inside the interior of the cell.
It calls coordinates-inside? which is defined for every cell class."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Reference cell dictionary
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; (display-ht *reference-cells*)
(defvar *reference-cells* (make-hash-table :test #'eq)
  "Table of all reference cells.")

(defun find-reference-cell-with (func)
  (let ((elems (find-reference-cells-with func)))
    (cond ((null elems) nil)
	  ((single? elems) (car elems))
	  (t (error "not unique")))))

(defun find-reference-cells-with (func)
  (loop for cell being each hash-key of *reference-cells*
	when (funcall func cell) collect cell))

(defun find-reference-cell-from-factors (factors)
  (find-reference-cell-with
   #'(lambda (cell)
       (equal (factor-simplices cell) factors))))

(defun find-reference-cell-from-factor-dimensions (factor-dims)
  (find-reference-cell-with
   #'(lambda (cell)
       (equal (mapcar #'dimension (factor-simplices cell))
	      factor-dims))))

(defun reference-cell? (cell)
  (and (gethash cell *reference-cells*) t))

(defun add-reference-cell! (refcell)
  (setf (gethash refcell *reference-cells*) t))

(defun remove-reference-cell! (refcell)
  (remhash refcell *reference-cells*))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Activation protocol of a cell-class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The activation of a cell-class is done in the following way:
;;; 1. Generate an cell-class containing the factor simplices.
;;; 2. Generate a reference cell for this class
;;; 3. Call initialize-cell-class on this reference cell.
;;;
;;; An exception is the initialization of the vertex class which is
;;; initialized by hand because it is needed by general refinement.

(defvar *cell-class-initialization-functions* ()
  "Register of cell-class initialization functions.")

(defun adjoin-cell-class-initialization-function (sym)
  "Add sym to the end of the initialization protocol, if it is not yet
included."
  (unless (member sym *cell-class-initialization-functions*)
    (push sym *cell-class-initialization-functions*)))

(defun setup-cell-class (refcell)
  "Initializes the cell-class slots from refcell and factor-simplices and
registers the reference element."
  (let ((cell-class (cell-class refcell)))
    (setf (cell-class-reference-cell cell-class) refcell)
    (setf (cell-class-dimension cell-class)
	  (if (zerop (length (boundary refcell)))
	      0
	      (1+ (dimension (aref (boundary refcell) 0)))))
    (setf (cell-class-nr-of-sides cell-class) (length (boundary refcell)))
    (create-subcell-info refcell)
    (setf (cell-class-nr-of-vertices cell-class) (length (subcells-of refcell 0)))))

(defun initialize-cell-class (refcell)
  "Initializes the cell-class by calling the initialization functions in
reverse order."
  (when (find-reference-cell-from-factors (factor-simplices refcell))
    (error "A reference cell with the same factors already exists."))
  (setup-cell-class refcell)
  ;; this is a potential problem
  (let ((initialization-done nil))
    (add-reference-cell! refcell)	; add temporarily
    (unwind-protect
	 (progn
	   (loop for sym in (reverse *cell-class-initialization-functions*)
		 do (funcall sym refcell))
	   (setq initialization-done t))
      (remove-reference-cell! refcell))
    (when initialization-done
      ;; add reference cell
      (add-reference-cell! refcell))))

;;; add this function to the initialization protocol
(adjoin-cell-class-initialization-function 'setup-cell-class)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Vertex class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The vertex class and reference class are generated here, because they
;;; are needed for the initialization protocol.  On the other hand, the
;;; initialization of the *cell-class* information is postponed so that
;;; actions on vertices are not yet possible.

(defstruct (<vertex> (:conc-name vertex-) (:include <cell>)
		     (:print-function print-vertex))
  "The basic (zero-dimensional) vertex.  This is a cell where the mapping
slot is usually filled with a position vector.  This class is used in
skeleton.lisp and is therefore defined immediately.")

(definline vertex-position (vtx) (mapping vtx))

#+(or)
(defclass <vertex> (<cell>)
  ()
  (:documentation "The basic (zero-dimensional) vertex.  This is a cell
where the mapping slot is usually filled with a position vector."))

(defvar *vertex-class* (make-<cell-class> :factor-simplices ())
  "The class of a vertex.")

(defun make-vertex (position)
  "The vertex constructor."
  (make-<vertex>
   :cell-class *vertex-class*
   :mapping (if (typep position 'double-vec)
		position
		(map 'double-vec #'identity position))))
  
(defvar *reference-vertex* (make-vertex (double-vec))
  "The reference vertex.")


;;; Testing: in mesh-tests.lisp

(defun test-cell ()
  "Tests are done when initializing the classes."
  (display-ht *reference-cells*)
  *cell-class-initialization-functions*
  )

(tests::adjoin-femlisp-test 'test-cell)
