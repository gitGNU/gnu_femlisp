;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; fe.lisp - Basic FE interface
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

(in-package :discretization)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <dof>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defstruct (<dof> (:conc-name dof-))
  "<dof> = degree of freedom.  It is defined as a functional defined by
integration over a sub-cell or by evaluation at a local coordinate of a
sub-cell of a reference cell.  <dof>-constructors are provided by
modules like lagrange.lisp, maybe later hermite.lisp.  The components
contain the following information:

index:           the index of the dof in the cell-dof array
subcell-index:   the index of the subcell on which the dof is defined
in-vblock-index: the index of the dof in the subcell vblock
subcell:         the subcell (information only)
coord:           the local coord of the dof (in the subcell)
gcoord:          the g-coord of the dof (on the cell, not the subcell)
functional:      an application to a function defined on the
                 (fe) reference-cell (not the subcell)
"
  (index -1 :type fixnum)
  (subcell-index -1 :type fixnum)
  (in-vblock-index -1 :type fixnum)
  (vblock-length -1 :type fixnum)
  (subcell (required-argument) :type <cell>)
  (coord (required-argument) :type double-vec)
  (gcoord (required-argument) :type double-vec)
  (functional nil))

(defmethod evaluate ((dof <dof>) fe-func)
  (evaluate (dof-functional dof) fe-func))

(defun interior-dof? (dof) (zerop (dof-subcell-index dof)))

(defstruct (<vector-dof> (:include <dof>) (:conc-name dof-))
  "A dof of a vector finite element.  Has an additional component field."
  (component 0 :type positive-fixnum))

(defun dof->vector-dof (dof component subcell-offsets)
  "Generates a vector-dof from a scalar dof."
  (make-<vector-dof>
   :index (dof-index dof)
   :subcell-index (dof-subcell-index dof)
   :in-vblock-index (+ (aref subcell-offsets (dof-subcell-index dof))
		       (dof-in-vblock-index dof))
   :subcell (dof-subcell dof)
   :coord (dof-coord dof)
   :gcoord (dof-gcoord dof)
   :functional (dof-functional dof)
   :component component))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <fe>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <fe> ()
  ((refcell :reader reference-cell :initarg :cell :type <cell>)
   (discretization :accessor discretization :initarg :discretization)
   (dofs :reader fe-dofs :initform () :initarg :dofs :type list)
   (basis :reader fe-basis :initform () :initarg :basis :type list)
   (properties :accessor properties :initform () :type list))
  (:documentation "A finite element <fe> is given for each reference cell,
e.g. <2-simplex>.  dofs are the degrees of freedom associated with the
cell, basis is the dual basis to dofs in some polynomial space.
subcell-ndofs is the number of ndofs on each subcell.  subcell-indices is a
list of indices for all subcells with dofs.  Usually, the <fe> will occur
as values of a procedure or as values in a hash-table with the reference
cells as keys."))

(defmethod nr-of-dofs ((fe <fe>)) (length (fe-basis fe)))
(defmethod nr-of-inner-dofs ((fe <fe>)) (aref (subcell-ndofs fe) 0))
(defmethod nr-of-components ((fe <fe>)) 1)
(defmethod components ((fe <fe>)) (vector fe))

(definline subcell-ndofs (fe)
  (the fixnum-vec (getf (properties fe) 'SUBCELL-NDOFS)))
(definline subcell-indices (fe)
  (the list (getf (properties fe) 'SUBCELL-INDICES)))
(definline inner-dof-indices (fe)
  (the list (getf (properties fe) 'INNER-DOF-INDICES)))

(defmacro do-dofs ((dof fe) &body body)
  `(loop for ,dof of-type <dof> in (fe-dofs ,fe) do
    ,@body))

(defmethod initialize-instance :after ((fe <fe>) &key &allow-other-keys)
  (with-slots (refcell properties)
    fe
    (let ((subcell-ndofs (make-fixnum-vec (nr-of-subcells refcell) 0)))
      (do-dofs (dof fe)
	(incf (aref subcell-ndofs (dof-subcell-index dof))))
      (setf (getf properties 'SUBCELL-NDOFS) subcell-ndofs)
      (setf (getf properties 'SUBCELL-INDICES)
	    (loop for nr-dofs across subcell-ndofs and i from 0
		  unless (zerop nr-dofs) collect i)))
    (setf (getf properties 'INNER-DOF-INDICES)
	  (coerce (loop for i below (nr-of-inner-dofs fe)
			and dof in (fe-dofs fe) collect
			(dof-in-vblock-index dof))
		  'vector))))

(defmethod make-local-vec ((fe <fe>) &optional (multiplicity 1))
  (make-real-matrix (nr-of-dofs fe) multiplicity))

(defmethod make-local-mat ((fe <fe>))
  (make-real-matrix (nr-of-dofs fe)))

(defmethod interpolate-on-refcell ((fe <fe>) function)
  "Interpolates FUNC on the reference cell of FE using FE."
  (let ((values (loop for dof in (subseq (fe-dofs fe) 0 (nr-of-inner-dofs fe))
		      collecting (evaluate dof function))))
    (when values
      (assert (apply #'= (mapcar #'multiplicity values)))
      (let ((vblock (make-real-matrix (nr-of-inner-dofs fe)
				      (multiplicity (first values)))))
	(loop for value in values and i from 0 do
	      (if (numberp value)
		  (setf (vref vblock i) value)
		  (minject vblock value i 0)))
	vblock))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <vector-fe>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <vector-fe> ()
  ((refcell :reader reference-cell :type <cell>)
   (discretization :accessor discretization :initarg :discretization)
   (components :accessor components :initarg :components :type simple-vector)
   (dofs :accessor fe-dofs :type list)
   (properties :reader properties :initform () :type list))
  (:documentation "Finite element for vector functions.  Components is a
vector of scalar finite elements.  Local-offset is an array of the same
length which contains the offsets for each component in the local
discretization vector.  Subcell-offsets is an array consisting of arrays
which yield such an offset for every subcell."))

(defmethod nr-of-dofs ((fe <vector-fe>))
  (reduce #'+ (components fe) :key #'nr-of-dofs))
(defmethod nr-of-inner-dofs ((fe <vector-fe>))
  (reduce #'+ (components fe) :key #'nr-of-inner-dofs))
(defmethod nr-of-components ((fe <vector-fe>))
  (length (components fe)))

;;; property content
(definline local-offset (fe)
  "Reader for the local-offset of this component in the local cell vector."
  (the fixnum-vec (getf (properties fe) 'LOCAL-OFFSET)))

(definline subcell-offsets (fe)
  "Reader for subcell-offsets.  This is an array of length the number of
components.  Each component is an array giving the offset of this component
in a sparse vector value block corresponding to the subcell."
  (the simple-vector (getf (properties fe) 'SUBCELL-OFFSETS)))

(defmethod initialize-instance :after ((vecfe <vector-fe>) &key &allow-other-keys)
  (with-slots (components dofs properties)
    vecfe
    (assert components)
    (assert (same-p components :key #'reference-cell))
    ;; setup derived slots
    (setf (slot-value vecfe 'refcell)
	  (reference-cell (aref components 0)))
    (let* ((refcell (reference-cell vecfe))
	   (nr-comps (length components))
	   (nr-subcells (nr-of-subcells refcell))
	   (local-offset (make-array nr-comps :element-type 'fixnum))
	   (subcell-offsets (make-array nr-comps))
	   (subcell-ndofs (make-array nr-subcells :element-type 'fixnum)))
      ;; fill arrays: local-offset, subcell-offsets, subcell-ndofs
      (loop with local-off = 0
	    and subcell-off = (make-fixnum-vec nr-subcells 0)
	    for i from 0
	    and fe across components do
	    (setf (aref local-offset i) local-off)
	    (setf (aref subcell-offsets i) (copy-seq subcell-off))
	    (incf local-off (nr-of-dofs fe))
	    (dotimes (j nr-subcells)
	      (incf (aref subcell-ndofs j)
		    (aref (subcell-ndofs fe) j))
	      (incf (aref subcell-off j)
		    (aref (subcell-ndofs fe) j))))
      (setf (getf properties 'LOCAL-OFFSET) local-offset)
      (setf (getf properties 'SUBCELL-OFFSETS) subcell-offsets)
      (setf (getf properties 'SUBCELL-NDOFS) subcell-ndofs)
      ;; set subcell-indices
      (setf (getf properties 'SUBCELL-INDICES)
	    (loop for nr-dofs across (subcell-ndofs vecfe)
		  and i from 0
		  unless (zerop nr-dofs) collect i))
      ;; setup dofs
      (setq dofs
	    (loop for fe across components and comp from 0
		  and local-off across local-offset
		  and subcell-offset across subcell-offsets nconcing
		  (loop for dof in (fe-dofs fe) collecting
			(dof->vector-dof dof comp subcell-offset))))
      )))

(defmethod make-local-vec ((vecfe <vector-fe>) &optional (multiplicity 1))
  (map 'simple-vector
       #'(lambda (fe) (make-real-matrix (nr-of-dofs fe) multiplicity))
       (components vecfe)))

(defmethod make-local-mat ((vecfe <vector-fe>))
  (let* ((n-comps (nr-of-components vecfe))
	 (result (make-array (list n-comps n-comps) :initial-element nil)))
    (loop for i from 0 and fe1 across (components vecfe) do
	  (loop for j from 0 and fe2 across (components vecfe) do
		(setf (aref result i j)
		      (make-real-matrix (nr-of-dofs fe1) (nr-of-dofs fe2)))))
    result))

(defmethod interpolate-on-refcell ((vecfe <vector-fe>) function)
  "Interpolates FUNC on the reference cell of VECTOR-FE using VECTOR-FE."
  (error "NYI.  Do you really need it?")
  #+(or)(coerce
   (loop for fe in (components vecfe) and comp from 0 collect
	 (interpolate-on-refcell
	  fe (sliced-function function comp)))
   'vector))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <fe-discretization>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <fe-discretization> (<discretization>)
  ()
  (:documentation "FE discretization base class."))

(defgeneric get-fe (fe-disc cell)
  (:documentation "Returns the finite element for the given discretization
and reference cell."))

(defgeneric quadrature-rule (fedisc fe)
  (:documentation "Computes the quadrature rule to be used for this finite element."))

(defmethod quadrature-rule ((fedisc <fe-discretization>) fe)
  "Standard quadrature rule for fe."
  (let ((refcell (reference-cell fe))
	(order (discretization-order fedisc)))
    (gauss-rule (mapcar #'dimension (factor-simplices refcell))
		(if (typep refcell '<simplex>)
		    order
		    (1+ order)))))

(defclass <standard-fe-discretization> (<fe-discretization>)
  ((cell->fe :initarg :cell->fe))
  (:documentation "For this class the finite elements are obtained from a
cell->fe mapping given as a class slot."))

(defmethod get-fe ((disc <standard-fe-discretization>) cell)
  (funcall (slot-value disc 'cell->fe)
	   (reference-cell cell)))

(defclass <scalar-fe-discretization> (<standard-fe-discretization>)
  ((order :reader discretization-order :initarg :order))
  (:documentation "Class for scalar fe discretizations."))

(defclass <vector-fe-discretization> (<standard-fe-discretization>)
  ((components :accessor components :initarg :components
	       :type (vector <scalar-fe-discretization>)))
  (:documentation "Vector FE discretization class."))

(defmethod initialize-instance ((disc <vector-fe-discretization>) &key components)
  "Combines scalar fe discretization to form a vector fe discretization."
  (setf (slot-value disc 'cell->fe)
	(compose
	 (memoize-1
	  #'(lambda (refcell)
	      (make-instance
	       '<vector-fe> :cell refcell :discretization disc :components
	       (map 'vector #'(lambda (comp-disc) (get-fe comp-disc refcell))
		    components))))
	   #'reference-cell))
  (setf (slot-value disc 'components)
	(coerce (loop for comp-disc in components and i from 0 collect
		      (make-instance
		       '<scalar-fe-discretization>
		       :order (discretization-order comp-disc)
		       :cell->fe
		       #'(lambda (cell)
			   (aref (components (get-fe disc cell)) i))))
		'vector))
    disc)


(defmethod discretization-order ((vecfe-disc <vector-fe-discretization>))
  (loop for disc across (components vecfe-disc)
	maximize (discretization-order disc)))

(defmethod nr-of-components ((vecfe-disc <vector-fe-discretization>))
  (length (components vecfe-disc)))

(defmethod component ((fedisc <standard-fe-discretization>) i)
  (make-instance
   '<standard-fe-discretization>
   :order (discretization-order fedisc)  ; might be wrong!
   :cell->fe
   #'(lambda (cell)
       (aref (components (funcall (slot-value fedisc 'cell->fe) cell))
	     i))))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ip-values, ip-gradients
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The following two functions might be the standard way to generate
;;; the quadrature information given an fe and a quadrature rule.
;;; They should be used in their memoized form

(defun base-function-values-at-ips (fe qrule)
  "Returns a list of nr-ip float-matrices of dimension (n-basis x 1)."
  (map 'vector
       #'(lambda (ip)
	   (make-real-matrix
	    (loop for shape in (fe-basis fe)
		  collect (list (evaluate shape (ip-coords ip))))))
        (integration-points qrule)))
(memoize-symbol 'base-function-values-at-ips :test 'equal)

(defun base-function-gradients-at-ips (fe qrule)
  "Returns a list of nr-ip float-matrices of dimension (n-basis x dim)."
  (map 'vector
       #'(lambda (ip)
	   (make-real-matrix
	    (loop for shape in (fe-basis fe) collect
		  (evaluate-gradient shape (ip-coords ip)))))
       (integration-points qrule)))
(memoize-symbol 'base-function-gradients-at-ips :test 'equal)

(defun vector-fe-ip-values (vecfe qrule)
  (apply #'map 'vector #'vector
	 (map 'list (rcurry 'ip-values qrule)
	      (components vecfe))))
(memoize-symbol 'vector-fe-ip-values :test 'equal)

(defun vector-fe-ip-gradients (vecfe qrule)
  (apply #'map 'vector #'vector
	 (map 'list (rcurry 'ip-gradients qrule)
	      (components vecfe))))
(memoize-symbol 'vector-fe-ip-gradients :test 'equal)

(defmethod ip-values ((fe <fe>) qrule)
  "Returns a list of nr-ip float-matrices of dimension (n-basis x 1)."
  (funcall #'base-function-values-at-ips fe qrule))

(defmethod ip-gradients ((fe <fe>) qrule)
  (funcall #'base-function-gradients-at-ips fe qrule))

(defmethod ip-values ((fe <vector-fe>) qrule)
  "Returns a list of vectors of length components."
  (funcall #'vector-fe-ip-values fe qrule))

(defmethod ip-gradients ((fe <vector-fe>) qrule)
  (funcall #'vector-fe-ip-gradients fe qrule))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Polynomial spaces on cells
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun P-nomials-of-degree (simplex deg &optional (type '=))
  "Returns a list of monomials of degree = deg or <= deg for a simplex."
  (n-variate-monomials-of-degree (dimension simplex) deg type))

(defun Q-nomials-of-degree (cell deg &optional (type '=))
  "Builds the Qn = Pn-Pn-Pn ... on a tensorial cell."
  (cond ((or (vertex-p cell) (simplex-p cell))
	 (P-nomials-of-degree cell deg type))
	((eq type '<=)
	 (apply #'map-product
		#'(lambda (&rest polys)
		    (make-polynomial
		     (reduce #'poly* polys)))
		;; later factors must be shifted
		(loop for factor in (factor-simplices cell)
		      for offset-dim = 0 then (+ offset-dim (dimension factor))
		      collecting
		      (loop for poly in (P-nomials-of-degree factor deg '<=)
			    collecting
			    (encapsulate (coefficients poly) offset-dim)))))
	(t (remove-if #'(lambda (poly) (< (maximal-partial-degree poly) deg))
		      (Q-nomials-of-degree cell deg '<=) ))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; gram matrix computation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun gram-matrix (vectors functionals pairing)
  "Computes the Gram matrix of vectors and functionals wrt pairing."
  (make-real-matrix
   (loop for phi in vectors
	 collect (loop for psi in functionals
		       collect (funcall pairing psi phi)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; dual basis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun compute-duals (vectors functionals pairing)
  "Combines a set of vectors such that they are dual to a set of
functionals wrt a given pairing."
  (unless (= (length vectors) (length functionals))
    (error "number of vectors and functionals must agree"))
  (let ((inverse-gram (m/ (gram-matrix vectors functionals pairing))))
    (loop for row from 0 below (nrows inverse-gram)
	  collect
	  (loop with sum = (zero (first vectors))
		for col from 0 below (ncols inverse-gram)
		and vec in vectors do
		(m+! (scal (mref inverse-gram row col) vec) sum)
		finally (return sum)))))

(defun projection-coefficients (vector basis pairing)
  "Returns the coefficients for a basis representation of a projection of vector
on the subspace given by basis wrt the scalar product pairing."
  (gesv! (gram-matrix basis basis pairing)
	 (make-real-matrix
	  (mapcar #'(lambda (basis-vector) (funcall pairing basis-vector vector))
		  basis))))

(defun project (vector basis pairing)
  "Projects vector on the subspace given by basis orthogonal to the
scalar product in pairing."
  (loop with result = vector
	with coeffs = (and basis (projection-coefficients vector basis pairing))
	for phi in basis
	and i from 0
	do (m-! (scal (vref coeffs i) phi) result)
	finally (return result)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; fe-cell-geometry
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun fe-cell-geometry (cell qrule)
  "Collects cell geometry information inside a property list."
  (loop for ip in (integration-points qrule)
	for local-coord = (ip-coords ip)
	for Dphi = (local->Dglobal cell local-coord)
	for volume = (area-of-span Dphi)
	collect local-coord into local-coords
	collect (local->global cell local-coord) into global-coords
	collect Dphi into gradients
	collect volume into volumes
	collect (when (= (nrows Dphi) (ncols Dphi)) (m/ Dphi)) into gradient-inverses
	collect (* (ip-weight ip) volume) into weights
	finally
	(return (list :cell cell :local-coords local-coords :global-coords global-coords
		      :gradients gradients :volume volumes
		      :gradient-inverses gradient-inverses :weights weights))))

;;; Testing: (test-fe)

(defun test-fe ()
  (Q-nomials-of-degree *unit-interval* 1 '=)
  (Q-nomials-of-degree *reference-vertex* 0 '<=)
  (Q-nomials-of-degree *unit-cube* 2)
  #+(or) ; works only if lagrange.lisp has been evaluated
  (let ((gm (gram-matrix (Q-nomials-of-degree *unit-interval* 1 '<=)
			 (lagrange-dofs *unit-interval* 1) #'evaluate)))
    (assert (< (abs (- 1.0 (mref (m* (m/ gm) gm) 0 0)) ) 1.0e-10)))
  )

(fl.tests:adjoin-test 'test-fe)
