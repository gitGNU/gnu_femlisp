;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; lagrange.lisp - Lagrange Finite Elements
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

;;;; This module provides Lagrange dofs and Lagrange basis functions
;;;; of arbitrary order and for arbitrary tensorial cells.

(defun lagrange-coords-1d (order type)
  (ecase type
    (:uniform
     (coerce (loop for i upto order
		collect (float (/ i order) 1.0)) 'vector))
    ;; The following choice of nodal points does not work for simplices with
    ;; dim>=3 and order>=3 because nodal points on the sides do not fit.  But
    ;; it is much better suited for interpolation in the case of tensorial
    ;; elements.
    (:gauss-lobatto
     (coerce (gauss-lobatto-points-on-unit-interval (1- order)) 'vector))))
  
(defmethod lagrange-inner-coords ((vtx <vertex>) order type)
  (list (double-vec)))

(defmethod lagrange-inner-coords ((simplex <simplex>) order type)
  (let ((coords-1d (lagrange-coords-1d order type))
	(dim (dimension simplex))
	(result ()))
    (multi-for (x (make-fixnum-vec dim 1) (make-fixnum-vec dim (1- order)))
      (when (< (reduce #'+ x) order)
	(push (map 'double-vec #'(lambda (i) (aref coords-1d i)) x)
	      result)))
    (reverse result)))

(defmethod lagrange-inner-coords ((cell <tensorial>) order type)
  (apply #'map-product #'(lambda (&rest args) (apply #'concatenate 'double-vec args))
	 (mapcar #'(lambda (simplex) (lagrange-inner-coords simplex order type))
		 (factor-simplices cell))))

(defun lagrange-dofs (cell order type)
  (let ((lagrange-coords (lagrange-coords-1d order type)))
    (loop with dof-index = -1
	  for subcell across (subcells cell)
	  and i from 0
	  nconcing
	  (loop with subcell-coords = (lagrange-inner-coords subcell order type)
		for local in subcell-coords
		and j from 0 collect
		;; we need below that the coords are eql to the lobatto
		;; coords without any rounding error
		(let ((g (map 'double-vec
			      #'(lambda (coord)
				  (find-if #'(lambda (coord2) (< (abs (- coord coord2)) 1.0e-10))
					   lagrange-coords))
			      (l2g subcell local))))
		  (make-<dof>
		   :index (incf dof-index)
		   :subcell subcell
		   :subcell-index i
		   :in-vblock-index j
		   :vblock-length (length subcell-coords)
		   :coord local
		   :gcoord g
		   :functional
		   #'(lambda (func) (evaluate func g))))))))
(memoize-symbol 'lagrange-dofs)

(defun lagrange-basis-simplex (cell order type)
  "Computes the Lagrange basis for a cell.  Should be applied only for
simplices, because for tensorials the basis can be computed as a tensor
product which is faster."
  (mapcar #'eliminate-small-coefficients
	  (compute-duals (Q-nomials-of-degree cell order '<=)
			 (lagrange-dofs cell order type)
			 #'evaluate)))

(defun shapes-and-dof-coords (factor-simplices order type)
  "Computes simpulataneously shapes and dof-coords for a tensorial as a
tensor product."
  (if (null factor-simplices)
      (values (list (double-vec)) (list (make-polynomial '(1.0))))
      (multiple-value-bind (coords shapes)
	  (shapes-and-dof-coords (cdr factor-simplices) order type)
	(let* ((factor (car factor-simplices))
	       (f-dim (dimension factor))
	       (f-shapes (lagrange-basis-simplex factor order type))
	       (f-dofs (lagrange-dofs factor order type))
	       (product
		(loop for coord in coords and shape in shapes nconcing
		      (loop for f-dof in f-dofs and f-shape in f-shapes
			    for f-coord = (dof-gcoord f-dof) collecting
			    (cons (concatenate 'double-vec f-coord coord)
				  (poly* (shift-polynomial shape f-dim) f-shape))))))
	  (values (mapcar #'car product) (mapcar #'cdr product))))))

(defun lagrange-basis (cell order type)
  "Computes the Lagrange basis for a tensorial accelerated.  The idea is to
construct the shapes with their associated dof-coordinates as a product of
lower-dimensional shapes and coordinates."
  (multiple-value-bind (coords shapes)
      (shapes-and-dof-coords (factor-simplices cell) order type)
    (let ((table (make-hash-table :test #'equalp)))
      (loop for coord in coords and shape in shapes do
	    (setf (gethash coord table) shape))
      (loop for dof in (lagrange-dofs cell order type)
	    for shape = (gethash (dof-gcoord dof) table)
	    do (assert shape)
	    collecting shape))))
(memoize-symbol 'lagrange-basis)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; fe-class definitions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <lagrange-fe> ()
  ((type :initarg :type))
  (:documentation "Mixin for finite elements of Lagrange type."))

(defclass <scalar-lagrange-fe> (<scalar-fe-discretization> <lagrange-fe>)
  ()
  (:documentation "Scalar Lagrange FE."))

(defclass <vector-lagrange-fe> (<vector-fe-discretization> <lagrange-fe>)
  ()
  (:documentation "Vector Lagrange FE."))

(defun lagrange-fe (order &key nr-comps integration-order (type :uniform))
  "Constructor for Lagrange fe."
  (unless integration-order (setq integration-order order))
  (if nr-comps
      (make-instance
       '<vector-lagrange-fe>
       :components
       (make-list nr-comps :initial-element
		  (lagrange-fe order :integration-order integration-order :type type)))
      (let ((disc (make-instance '<scalar-lagrange-fe> :order order :type type)))
	(setf (slot-value disc 'cell->fe)
	      (compose
	       (memoize-1
		#'(lambda (refcell)
		    (assert (reference-cell-p refcell))
		    (make-instance
		     '<scalar-fe> :cell refcell :discretization disc
		     :dofs (lagrange-dofs refcell order type)
		     :basis (lagrange-basis refcell order type)
		     :order order)))
	       #'reference-cell))
	disc)))
(memoize-symbol 'lagrange-fe)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Isoparametric stuff
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun lagrange-basis-boundary (cell order type)
  (loop for phi in (lagrange-basis cell order type)
	for dof in (lagrange-dofs cell order type)
	unless (interior-dof? dof) collect phi))

(defun lagrange-basis-inner (cell order type)
  (loop for phi in (lagrange-basis cell order type)
	for dof in (lagrange-dofs cell order type)
	when (interior-dof? dof) collect phi))

(defun lagrange-boundary-dofs (cell order type)
  (loop	for dof in (lagrange-dofs cell order type)
	unless (interior-dof? dof) collect dof))

(defun lagrange-reference-parameters (refcell order type)
  "Computes an energy-minimizing extension to the interior for the
boundary lagrangian."
  (assert (reference-cell-p refcell))
  (let* ((fe-class (lagrange-fe order :type type))
	 (fe (get-fe fe-class refcell))
	 (nr-inner-dofs (nr-of-inner-dofs fe))
	 (result (make-hash-table)))
    (cond
      ((zerop nr-inner-dofs)
       (loop+ ((dof (fe-dofs fe)) (shape (fe-basis fe)))
	  do (setf (gethash dof result) shape)))
      (t
       (let* ((nr-dofs (nr-of-dofs fe))
	      (energy-mat (make-real-matrix nr-dofs))
	      (inner-indices (range< 0 nr-inner-dofs))
	      (boundary-indices (range< nr-inner-dofs nr-dofs)))
	 (loop with qrule = (quadrature-rule fe)
	       for ip in (integration-points qrule)
	       for gradients across (ip-gradients fe qrule) do
	       (gemm! (ip-weight ip) gradients gradients 1.0 energy-mat :nt))
	 (let* ((A_II (submatrix energy-mat :row-indices inner-indices :col-indices inner-indices))
		(A_IB (submatrix energy-mat :row-indices inner-indices :col-indices boundary-indices))
		(corr-mat (m* (m/ A_II) A_IB))
		(inner-shapes (lagrange-basis-inner refcell order type)))
	   (loop for j from 0
		 and dof in (lagrange-boundary-dofs refcell order type)
		 and phi in (lagrange-basis-boundary refcell order type) do
		 (setf (gethash dof result)
		       (loop for i from 0 and psi in inner-shapes do
			     (setf phi (axpy (- (mref corr-mat i j)) psi phi))
			     finally (return phi))))))))
    result))

(defun lagrange-polynomial-vector (cell order type)
  "Returns a vector of polynomials representing the isoparametric mapping."
  (loop with result = (copy (make-array (embedded-dimension cell)
					:initial-element (make-polynomial '(0))))
	with basis = (lagrange-reference-parameters (reference-cell cell) order type)
	with subcells = (subcells cell)
	for dof being the hash-keys of basis
	and phi being the hash-values of basis do
	(loop for comp across (local->global (aref subcells (dof-subcell-index dof))
					     (dof-coord dof))
	      and i from 0 do
	      (axpy! comp phi (aref result i)))
	finally (return result)))

(defun lagrange-mapping (order &optional (type :uniform))
  "Returns a function which maps a cell by a polynomial which is obtained
by interpolating the boundary map via Lagrange interpolation."
  #'(lambda (cell)
      (let ((poly-vec (lagrange-polynomial-vector cell order type)))
	(dbg :lagrange "Lagrange-map(~A)=~%~A" cell poly-vec)
	(make-instance
	 '<special-function>
	 :domain-dimension (dimension cell)
	 :image-dimension (length poly-vec)
	 :evaluator
	 #'(lambda (lcoord)
	     (loop with result = (make-double-vec (length poly-vec))
		   for poly across poly-vec
		   and i from 0 do
		   (setf (aref result i) (float (evaluate poly lcoord) 1.0))
		   finally (return result)))
	 :gradient
	 #'(lambda (lcoord)
	     (make-real-matrix
	      (loop for poly across poly-vec collect
		    (loop for i from 0 below (length lcoord) collect
			  (evaluate (differentiate poly i) lcoord)))))))))



;;; Testing
(defun test-lagrange ()
  (flet ((check-fe (fe-class refcell)
	   (let* ((fe (get-fe fe-class refcell))
		  (qrule (quadrature-rule fe)))
	     (integration-points qrule)
	     (list (ip-values fe qrule) (ip-gradients fe qrule)))))
    (check-fe (lagrange-fe 1) *reference-vertex*)
    (check-fe (lagrange-fe 3 :nr-comps 2) *unit-cube*))
  (lagrange-basis *reference-vertex* 1 :uniform)
  (lagrange-basis *unit-interval* 1 :uniform)
  (lagrange-basis *unit-quadrangle* 1 :uniform)
  (lagrange-dofs *unit-interval* 1 :uniform)
  (lagrange-dofs *unit-interval* 2 :uniform)
  (lagrange-basis *unit-triangle* 2 :uniform)
  (lagrange-dofs *unit-triangle* 2 :uniform)
  (lagrange-dofs *unit-triangle* 3 :uniform)
  (lagrange-dofs *unit-tetrahedron* 3 :uniform)
  (lagrange-basis *unit-triangle* 1 :uniform)
  (lagrange-dofs *unit-cube* 2 :uniform)
  (time (let ()
	  (lagrange-basis-simplex *unit-triangle* 15 :uniform)
	  nil))  ; 1.0-1.2s (toba)
  (lagrange-inner-coords *unit-quadrangle* 3 :uniform)
  (time (lagrange-reference-parameters *unit-quadrangle* 7 :uniform))

  ;; Isoparametric stuff
  (let ((center (make-vertex #(0.0 -4.0)))
	(east-vtx (make-vertex #(-1.0 1.0)))
	(north-vtx (make-vertex #(0.0 1.0))))
    ;; line segments
    (let ((seg-ce (make-line center east-vtx))
	  (seg-cn (make-line center north-vtx))
	  (seg-en (make-line east-vtx north-vtx)))
      (let ((tri (make-simplex (vector seg-en seg-cn seg-ce)))
	    (lcoord #(0.2 0.3)))
	(princ (det (local->Dglobal tri lcoord))) (terpri)
	(princ (local->global tri lcoord))
	(evaluate (funcall (lagrange-mapping 4) tri) lcoord))))

  (assert (eq (reference-cell (get-fe (lagrange-fe 1) *unit-interval*))
	      *unit-interval*))
  (setq *print-matrix* 7)
  (let ((fe (get-fe (lagrange-fe 1) *unit-interval*)))
    (fe-basis fe))
  (get-fe (lagrange-fe 3) *unit-interval*)
  (get-fe (lagrange-fe 3) (n-simplex 3))

  ;;
  (describe
   (make-instance
    '<vector-fe>
    :components (make-array 2 :initial-element (get-fe (lagrange-fe 1) *unit-interval*))))
  ;;
  (dbg-on :lagrange)
  (let* ((domain (n-ball-domain 2))
	 (mesh (make-mesh-from-domain domain :parametric (lagrange-mapping 1)))
	 (refined (refine mesh)))
    (check refined))
  (dbg-off :lagrange)
  )

(fl.tests:adjoin-test 'test-lagrange)

;;; (fl.discretization::test-lagrange)

