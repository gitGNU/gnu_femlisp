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

(in-package :fl.mesh)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Cell class definition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <cell> ()
  ()
  (:documentation "The basic cell class."))

(deftype cell-vec () '(simple-array <cell> (*)))

(defclass <standard-cell> (<cell>)
  ((boundary :reader boundary :initarg :boundary :type cell-vec :documentation
	     "A vector of boundary cells."))
  (:documentation "The standard cell in Femlisp is defined via its boundary."))

(defgeneric local->global (cell local-pos)
  (:documentation "local->global checks if a mapping is given for the cell.
If yes, then this mapping is evaluated.  If no, then the function l2g is called
which should do a multilinear interpolation from the cell's corners."))

(defgeneric local->Dglobal (cell local-pos)
  (:documentation "local->Dglobal checks if a mapping is given for the cell.
If yes, then the gradient of this mapping is evaluated (if available).  If no,
then the function l2Dg is called which gives the gradient for a multilinear
interpolation from the cell's corners."))

(defgeneric midpoint (cell)
  (:documentation "Returns cell midpoint in global coordinates.")
  (:method (cell) "Default method"
	   (local->global cell (local-coordinates-of-midpoint cell))))

(defgeneric origin (cell)
  (:documentation "Returns cell origin in global coordinates.")
  (:method (cell) "Default method"
	   (local->global cell (make-double-vec (dimension cell)))))

(defparameter *print-cell* :midpoint
  "If set to :MIDPOINT, prints the midpoint of a cell.")

(defmethod print-object :after ((cell <cell>) stream)
  (case *print-cell*
    (:midpoint (format stream "{MP=~A}" (midpoint cell)))))

;;; Specializing for vertices is needed very soon.  Therefore, we specify
;;; the vertex class without initializing its class information
;;; immediately.
(defclass <vertex> (<cell>)
  ((position :reader vertex-position :initarg :position :type double-vec))
  (:documentation "The vertex class."))

(declaim (inline vertex? vertex-p))
(defun vertex? (cell) (typep cell '<vertex>))
(defun vertex-p (cell) (typep cell '<vertex>))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; class information for cells
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (cell-class-information (:conc-name ci-))
  "Structure containing information necessary for cell handling."
  (dimension -1 :type fixnum)
  (reference-cell nil)
  (factor-simplices ())
  (nr-of-sides -1 :type fixnum)
  (nr-of-vertices -1 :type fixnum)
  (nr-of-subcells -1 :type fixnum)
  (boundary-indices-of-subcells nil)
  (refine-info nil)
  (refinement-rules #()))

;;; for debugging purposes this is preferable
(declaim (notinline cell-class-information
		 cell-dimension reference-cell factor-simplices
		 nr-of-sides nr-of-vertices nr-of-subcells
		 boundary-indices-of-subcells refine-info))

(defun cell-class-information (obj)
  "Returns the cell information for the class which is stored as a property
of the class symbol."
  (get (etypecase obj
	 (symbol obj)
	 (<cell> (class-name (class-of obj)))
	 (standard-class (class-name obj)))
       'cell-class-information))

(defun (setf cell-class-information) (value obj)
  "This function should only be used during class initialization."
  (setf (get (etypecase obj
	       (symbol obj)
	       (<cell> (class-name (class-of obj)))
	       (standard-class (class-name obj)))
	     'cell-class-information)
	value))

(defmacro with-ci-slots (slots ci &body body)
  "Multiple and write access to cell-information slots."
  `(symbol-macrolet
    ,(loop for item in slots
	   for accessor = `(,(intern (concatenate 'string "CI-" (symbol-name item))) ,ci)
	   collect `(,item ,accessor))
    ,@body))

(defmacro with-cell-class-information (slots cell-class &body body)
  "Multiple and write access to cell-information slots."
  (with-gensyms (ci)
    `(let ((,ci (cell-class-information ,cell-class)))
      (with-ci-slots ,slots ,ci ,@body))))

(defmacro with-cell-information (slots cell &body body)
  (with-gensyms (class)
    `(let ((,class (class-of ,cell)))
      (with-cell-class-information ,slots ,class ,@body))))

(defun factor-simplices (cell)
  "Returns the factor-simplices."
  (with-cell-information (factor-simplices) cell factor-simplices))
(defun nr-of-sides (cell)
  "Returns the number of boundary faces."
  (with-cell-information (nr-of-sides) cell nr-of-sides))
(defun nr-of-vertices (cell)
  "Returns the number of vertices."
  (with-cell-information (nr-of-vertices) cell nr-of-vertices))
(defun nr-of-subcells (cell)
  "Returns the number of subcells."
  (with-cell-information (nr-of-subcells) cell nr-of-subcells))
(defun refine-info (cell)
  "Returns refinement information for the cell."
  (with-cell-information (refine-info) cell refine-info))

(defun refinement-rules (cell)
  "Returns refinement information for the cell."
  (with-cell-information (refinement-rules) cell refinement-rules))

;;; because of their nice names, the following are implemented as generic
;;; functions
(defmethod dimension ((cell <cell>))
  "Returns the dimension of the cell."
  (with-cell-information (dimension) cell dimension))

(defmethod reference-cell ((cell <cell>))
  "Returns the cell's or cell-classes reference-cell."
  (with-cell-information (reference-cell) cell reference-cell))

(defmethod reference-cell ((class standard-class))
  "Returns the cell information also when called for a cell class."
  (with-cell-class-information (reference-cell) class reference-cell))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Parametrized cell classes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <mapped-cell> ()
  ((mapping :accessor mapping :initarg :mapping))
  (:documentation "A mixin which distinguishes cells which are mapped by a
special mapping."))

(defclass <distorted-cell> () ()
  (:documentation "A mixin which distinguishes if the cell mapping is a
distortion of the multilinear mapping."))

(definline mapped-p (cell) (subtypep (class-of cell) '<mapped-cell>))

(defun mapped-cell-class (class &optional distorted)
  "Constructs a cell class with <mapped-cell> mixin."
  (assert (eq (symbol-package (if (symbolp class) class (class-name class)))
	      (find-package "FL.MESH")))
  (if (subtypep class '<mapped-cell>)
      class
      (let* ((unmapped-class (class-name class))
	     (mapped-class
	      (intern (concatenate 'string "<MAPPED-"
				   (subseq (symbol-name unmapped-class) 1))
		      "FL.MESH")))
	(or (find-class mapped-class nil)
	    (let ((new-class (eval `(defclass ,mapped-class
				     (,@(when distorted 'fl.mesh::<distorted-cell>)
					fl.mesh::<mapped-cell>
					,unmapped-class) ()))))
	      (setf (cell-class-information new-class)
		    (cell-class-information class))
	      new-class)))))

(defun unmapped-cell-class (class)
  "Returns the cell class without mapping mixin."
  (assert (eq (symbol-package (if (symbolp class) class (class-name class)))
	      (find-package "FL.MESH")))
  (if (subtypep class '<mapped-cell>)
      (let* ((mapped-class (class-name class))
	     (unmapped-class
	      (intern (concatenate 'string "<" (subseq (symbol-name mapped-class) 8))
		      "FL.MESH")))
	(find-class unmapped-class nil))
      class))

(defmethod mapping ((cell <cell>))
  "No mapping for ordinary cells."
  nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Further routines
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric vertices (cell)
  (:documentation "Returns a list of all vertices of the cell."))

(defgeneric copy-cell (cell)
  (:documentation "Copy constructor for cells.  It is guaranteed that the
cell boundary is copied."))

(definline reference-cell-p (cell)
  "Tests if a cell is a reference cell."
  (eq (reference-cell cell) cell))

(defmethod copy-cell ((cell <standard-cell>))
  "Copy constructor for cells."
  (make-instance (class-of cell) :boundary (copy-seq (boundary cell))))

(defmethod copy-cell ((cell <mapped-cell>))
  (let ((copy (call-next-method)))
    (setf (slot-value copy 'mapping)
	  (slot-value cell 'mapping))
    copy))

(defmethod check ((cell <cell>))
  (loop with side-dim = (1- (dimension cell))
	for side across (boundary cell)	do
	(assert (= side-dim (dimension side)))))

(defmethod check :after ((cell <mapped-cell>))
  (with-slots (mapping) cell
    ;; check if mapping is reasonable
    (assert (differentiable-p mapping))
    (let* ((mp (local-coordinates-of-midpoint cell))
	   (grad (evaluate-gradient mapping mp))
	   (numgrad (evaluate (numerical-gradient mapping) mp)))
      (assert (< (norm (m- grad numgrad))
		 (* 1.0e-2 (max (norm grad) (norm numgrad))))))))

(defmethod describe-all ((cell <cell>))
  (describe cell)
  (mapc #'describe-all (boundary cell)))

;;; Ordered access (i.e. the ordering depends only on the cell's class) to
;;; sub-cells

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Access to subcells (multiple cells are allowed for generality)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun subcell-access-indices (cells)
  "Constructs an ordering of all subcells of the list @arg{cells}.  All
cells in this list should have the same dimension."
  (let ((remaining (list->queue cells))
	(parent-indices (make-list (length cells) :initial-element -1))
	(boundary-indices (make-list (length cells) :initial-element -1)))
    (loop for parent-index from 0
	 for cell = (dequeue remaining) while cell do
	 (loop for side across (boundary cell) and k from 0
	      unless (member side (queue->list remaining)) do
	      (push parent-index parent-indices)
	      (push k boundary-indices)
	      (enqueue side remaining)))
    (values (reverse parent-indices) (reverse boundary-indices))))

(defun generate-subcells-method-source (refcell)
  "Generates code for a @function{subcells} method for @arg{cell-class}."
  (multiple-value-bind (parent-indices boundary-indices)
      (subcell-access-indices (list refcell))
    `(defmethod subcells ((cell ,(class-name (class-of refcell))))
       (declare (optimize speed))
       (let ((result (make-array ,(length parent-indices) :initial-element nil)))
	 ,@(loop with i = 0 until (= i (length parent-indices)) collect
		(let ((parent-index (elt parent-indices i))
		      (boundary-index (elt boundary-indices i)))
		  (if (= boundary-index -1)
		      (prog1
			  `(setf (aref result ,i) cell)
			(incf i))
		      `(let* ((cell (aref result ,parent-index))
			      (boundary (boundary cell)))
			 (declare (type cell-vec boundary))
			 ,@(loop for j from i below (length parent-indices)
			      for parent-index-2 = (elt parent-indices j)
			      while (= parent-index-2 parent-index) collect
			      `(setf (aref result ,j) (aref boundary ,(elt boundary-indices j)))
			      do (incf i))))))
	 result))))

(defgeneric subcells (cell)
  (:documentation
   "Returns a vector containing all subcells of a given cell.  The code is
special to each class and often automatically generated by
@function{generate-subcell-access-code}."))

(defun generate-subcells-method (refcell)
  "Defines a method for @arg{subcells}.  This function should be called the
cell class is initialized."
  (assert (reference-cell-p refcell))
  (fl.amop:compile-and-eval (generate-subcells-method-source refcell)))

(let (previous-cell previous-corners)
  (defmethod corners ((cell <cell>))
    "Returns the corners of the cell.  The last value is cached because it
is called often repeatedly on the same cell in l2g."
    (cond
      ((eq previous-cell cell) previous-corners)
      (t (setq previous-cell cell)
	 (setq previous-corners
	       (mapcar #'(lambda (vtx) (local->global vtx (double-vec)))
		       (vertices cell)))))))

(let (previous-cell previous-matrix)
  (defun corner-matrix (cell)
    "Returns a matrix whose columns are the corners of the cell."
    (cond
      ((eq previous-cell cell) previous-matrix)
      (t (setq previous-cell cell)
	 (setq previous-matrix
	       (let ((corners (corners cell)))
		 (make-instance
		  (standard-matrix 'double-float)
		  :nrows (length (car corners)) :ncols (length corners)
		  :store (apply #'concatenate 'double-vec corners))))))))

(defmethod g-corners ((cell <cell>))
  (mapcar #'(lambda (vtx) (l2g vtx (double-vec)))
	  (vertices cell)))

(defgeneric embedded-dimension (object)
  (:documentation "Dimension of the embedding space for object.")
  (:method ((cell <cell>))
    "Recursive definition, anchored at the definition for vertices."
    (embedded-dimension (aref (boundary cell) 0))))

(defmethod cell-mapping ((cell <cell>))
  "For non-mapped cells, this method returns a <special-function> which can
be called instead of @function{l2g} and @function{l2Dg}."
  ;; is the mapping linear?
  (let* ((local (make-double-vec (dimension cell)))
	 (A (l2Dg cell local)) (b (l2g cell local))
	 (linmap (make-instance '<linear-function> :A A :b b))
	 (delta (* (norm A) 1000 double-float-epsilon)))
    (if (every (lambda (corner refcell-corner)
		 (mzerop (m- corner (evaluate linmap refcell-corner)) delta))
	       (corners cell) (corners (reference-cell cell)))
	linmap
	(make-instance
	 '<special-function>
	 :domain-dimension (dimension cell) :image-dimension (embedded-dimension cell)
	 :evaluator (curry #'l2g cell) :gradient (curry #'l2Dg cell)))))

(defmethod cell-mapping ((cell <mapped-cell>))
  "Return the mappingFor non-mapped cells, this returns a <special-function> which is
equivalent to calling l2g and l2Dg."
  (slot-value cell 'mapping))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; local->global for isoparametric cells (as an after method)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod local->global ((cell <cell>) local-pos)
  (l2g cell local-pos))
(defmethod local->global ((cell <mapped-cell>) local-pos)
  (evaluate (slot-value cell 'mapping) local-pos))
(defmethod local->global ((cell <distorted-cell>) local-pos)
  (evaluate (slot-value cell 'mapping) (l2g cell local-pos)))

(defmethod local->Dglobal ((cell <cell>) local-pos)
  (l2Dg cell local-pos))
(defmethod local->Dglobal ((cell <mapped-cell>) local-pos)
  (evaluate-gradient (slot-value cell 'mapping) local-pos))
(defmethod local->Dglobal ((cell <distorted-cell>) local-pos)
  (let ((g (l2g cell local-pos))
	(Dg (l2Dg cell local-pos)))
    (m* (evaluate-gradient (slot-value cell 'mapping) g) Dg)))

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

(defvar *g2l-newton-steps* 20
  "Limit for the number of Newton steps for global->local.  Small numbers
should work for reasonable geometries and initial values.")

(defun simple-newton (f x-start &key approximate-gradient)
  "A simple Newton iteration for converting global to local coordinates."
  (loop with initial-defnorm = nil
     for defnorm = nil then new-defnorm
     for x = x-start then (m- x (gesv (funcall approximate-gradient x) defect))
     for defect = (funcall f x)
     for new-defnorm = (norm defect)
     for i = 0 then (1+ i)
     do
       (ensure initial-defnorm defnorm)
       (dbg :mesh "Step ~3D: |f| = ~8,3,2E, x = ~A~%" i new-defnorm x)
       (when (and defnorm (or (<= defnorm new-defnorm)
			      (<= defnorm (* 1.0e-10 initial-defnorm))))
	 (return x))
       (when (> i *g2l-newton-steps*)
	 (error "Newton did not converge fast enough.  Improve the mesh or
increase *g2l-newton-steps*."))))

(defmethod global->local ((cell <cell>) global-pos)
  "Does a Newton iteration or a Gauss-Newton method for approximating
global-pos by the cell mapping."
  (if (= (dimension cell) (length global-pos))
      (simple-newton #'(lambda (x) (m- (local->global cell x) global-pos))
		     (make-double-vec (dimension cell))
		     :approximate-gradient #'(lambda (x) (local->Dglobal cell x)))
      (let* ((A (local->Dglobal cell (make-double-vec (dimension cell))))
	     (At (transpose A))
	     (At*A (m* At A)))
	(simple-newton #'(lambda (x) (m* At (m- (local->global cell x) global-pos)))
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; global->embedded-local
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric global->embedded-local (cell global-pos)
  (:documentation "This computes a local coordinate which solves the
Ausgleichsproblem of mapping to a point as near as possible to global-pos.  It
involves more computational work than global->local.  As a second value, the
distance to global-pos is returned."))

(defmethod global->embedded-local ((cell <cell>) global-pos)
  (let ((x (simple-newton
	    #'(lambda (x)
		(m* (transpose (local->Dglobal cell x))
		    (m- (local->global cell x) global-pos)))
	    (make-double-vec (dimension cell))
	    :approximate-gradient
	    #'(lambda (x)
		(let ((A (local->Dglobal cell x)))
		  (m* (transpose A) A))))))
    (values x (norm (m- (local->global cell x) global-pos)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Initialization of a cell class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric generate-refine-info (refcell)
  (:documentation "Generates the refinement information for the reference
cell."))

(defgeneric skeleton (cell-or-cells)
  (:documentation "Returns a skeleton for the given cell or the given
cells."))

(defun initialize-cell-class (refcell factors)
  "Initializes the per-class information for the given reference cell."
  (let ((cell-class-information
	 (or (cell-class-information (class-of refcell))
	     (setf (cell-class-information (class-of refcell))
		   (make-cell-class-information
		    :reference-cell refcell :factor-simplices factors)))))
    (with-ci-slots (dimension reference-cell factor-simplices
			      nr-of-sides nr-of-subcells nr-of-vertices
			      boundary-indices-of-subcells
			      refinement-rules)
      cell-class-information
      (setq reference-cell refcell)
      (setq factor-simplices factors)
      (setq dimension (if (zerop (length (boundary refcell)))
			  0
			  (1+ (dimension (aref (boundary refcell) 0)))))
      (setq nr-of-sides (length (boundary refcell)))
      (generate-subcells-method refcell)
      (let ((refcell-subcells (subcells refcell)))
	(setf nr-of-subcells (length refcell-subcells))
	(setf boundary-indices-of-subcells
	      (map 'vector
		   #'(lambda (subcell)
		       (map 'vector (rcurry #'position refcell-subcells)
			    (boundary subcell)))
		   refcell-subcells))
	(setq nr-of-vertices
	      (count-if #'zerop refcell-subcells :key #'dimension)))
      ;; set the refinement info
      (generate-refine-info refcell)
      nil
      )))

;;; Testing: in mesh-tests.lisp

(defun test-cell ()
  (flet ((simple-linear-f (x) (+ 1.0 (* 0.5 x))))
    (simple-newton #'simple-linear-f 1.0 :approximate-gradient (constantly 0.5)))
  )

(fl.tests:adjoin-test 'test-cell)
