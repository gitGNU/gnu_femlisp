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

(in-package :fl.discretization)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; dof
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass dof ()
  ((index :reader dof-index :initarg :index :type fixnum
	  :documentation "index of the dof in the cell-dof array")
   (subcell-index :reader dof-subcell-index :initarg :subcell-index :type fixnum
		  :documentation "index of the reference subcell on which the dof is defined")
   (in-vblock-index :reader dof-in-vblock-index :initarg :in-vblock-index :type fixnum
		    :documentation "index of the dof in the subcell vblock")
   (subcell :reader dof-subcell :initarg :subcell
		  :documentation "reference subcell on which the dof is defined")
   (coord :reader dof-coord :initarg :coord :type double-vec
	    :documentation "local coordinate of the dof in the reference subcell")
   (gcoord :reader dof-gcoord :initarg :gcoord :type double-vec
	  :documentation "global coordinate of the dof on the reference cell")
   (functional :reader dof-functional :initarg :functional
	   :documentation "a functional for functions defined on the reference cell"))
  (:documentation
   "Degree of freedom in a finite element.  It is defined as a functional
defined by integration over a sub-cell or by evaluation at a local
coordinate of a sub-cell of a reference cell. "))

(defmethod evaluate ((dof dof) fe-func)
  (evaluate (dof-functional dof) fe-func))

(defun interior-dof? (dof) (zerop (dof-subcell-index dof)))

(defclass vector-dof (dof)
  ((component :reader dof-component :initarg :component :documentation
    "The component in the solution vector to which this @class{dof} belongs."))
  (:documentation "A dof of a vector finite element."))

(defun new-vector-dof-from-dof (dof component subcell-offsets)
  "Generates a vector-dof from a scalar dof."
  (make-instance 'vector-dof
		 :index (dof-index dof)
		 :subcell-index (dof-subcell-index dof)
		 :in-vblock-index (+ (aref subcell-offsets (dof-subcell-index dof))
				     (dof-in-vblock-index dof))
		 :subcell (dof-subcell dof)
		 :coord (dof-coord dof) :gcoord (dof-gcoord dof)
		 :functional (dof-functional dof)
		 :component component))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <fe>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <fe> (property-mixin)
  ((refcell :reader reference-cell :initarg :cell :type <cell>)
   (discretization :accessor discretization :initarg :discretization))
  (:documentation "Abstract base class for finite elements."))

(defclass <scalar-fe> (<fe>)
  ((dofs :reader fe-dofs :initform () :initarg :dofs :type list)
   (basis :reader fe-basis :initform () :initarg :basis :type list)
   (order :reader discretization-order :initarg :order))
  (:documentation "A finite element <fe> is given for each reference cell,
e.g. <2-simplex>.  dofs are the degrees of freedom associated with the
cell, basis is the dual basis to dofs in some polynomial space.
subcell-ndofs is the number of ndofs on each subcell.  subcell-indices is a
list of indices for all subcells with dofs.  Usually, the <scalar-fe> will occur
as values of a procedure or as values in a hash-table with the reference
cells as keys."))

(defmethod nr-of-dofs ((fe <scalar-fe>)) (length (fe-basis fe)))
(defmethod nr-of-inner-dofs ((fe <scalar-fe>)) (aref (subcell-ndofs fe) 0))
(defmethod nr-of-components ((fe <scalar-fe>)) 1)
(defmethod components ((fe <scalar-fe>)) (vector fe))

(defun subcell-ndofs (fe)
  (getf (properties fe) 'SUBCELL-NDOFS))
(defun subcell-indices (fe)
  (getf (properties fe) 'SUBCELL-INDICES))
(defun inner-dof-indices (fe)
  (getf (properties fe) 'INNER-DOF-INDICES))

(defmacro do-dof ((dof-and-shape fe &key (type :dof)) &body body)
  (with-gensyms (fe2)
    (destructuring-bind (dof shape)
      (ecase type
	(:dof (list dof-and-shape nil))
	(:shape (list nil dof-and-shape))
	(:dof-and-shape dof-and-shape))
      `(let ((,fe2 ,fe))
	 (map nil
	      (lambda ,(append (mklist dof) (mklist shape))
		,@body)
	      ,@(when dof `((fe-dofs ,fe2)))
	      ,@(when shape `((fe-basis ,fe2))))))))


(defmethod initialize-instance :after ((fe <scalar-fe>) &key &allow-other-keys)
  (with-slots (refcell properties)
    fe
    (let ((subcell-ndofs (make-fixnum-vec (nr-of-subcells refcell) 0)))
      (do-dof (dof fe)
	(incf (aref subcell-ndofs (dof-subcell-index dof))))
      (setf (getf properties 'SUBCELL-NDOFS) subcell-ndofs)
      (setf (getf properties 'SUBCELL-INDICES)
	    (loop for nr-dofs across subcell-ndofs and i from 0
		  unless (zerop nr-dofs) collect i)))
    (setf (getf properties 'INNER-DOF-INDICES)
	  (coerce (loop+ ((i (range :below (nr-of-inner-dofs fe)))
			  (dof (fe-dofs fe)))
		     collecting (dof-in-vblock-index dof))
		  'vector))))

(defmethod make-local-vec ((fe <scalar-fe>) &optional (multiplicity 1))
  (make-real-matrix (nr-of-dofs fe) multiplicity))

(defmethod make-local-mat ((fe <scalar-fe>))
  (make-real-matrix (nr-of-dofs fe)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; <vector-fe>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <vector-fe> (<fe>)
  ((components :accessor components :initarg :components :type simple-vector)
   (dofs :accessor fe-dofs :type list))
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

(defmethod component ((vecfe <vector-fe>) comp)
  (aref (components vecfe) comp))

(defmethod discretization-order ((vecfe <vector-fe>))
  (loop for fe across (components vecfe)
	maximize (discretization-order fe)))


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
  (declare (optimize safety debug))
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
	   (local-offset (make-fixnum-vec nr-comps))
	   (subcell-offsets (make-array nr-comps :initial-element nil))
	   (subcell-ndofs (make-fixnum-vec nr-subcells)))
      ;; fill arrays: local-offset, subcell-offsets, subcell-ndofs
      (loop with local-off = 0
	 and subcell-off = (make-fixnum-vec nr-subcells 0)
	 for i from 0 and fe across components
	 do
	   (setf (aref local-offset i) local-off)
	   (setf (aref subcell-offsets i) (copy-seq subcell-off))
	   (incf local-off (nr-of-dofs fe))
	   (m+! (subcell-ndofs fe) subcell-ndofs)
	   (m+! (subcell-ndofs fe) subcell-off))
      (setf (getf properties 'LOCAL-OFFSET) local-offset)
      (setf (getf properties 'SUBCELL-OFFSETS) subcell-offsets)
      (setf (getf properties 'SUBCELL-NDOFS) subcell-ndofs)
      ;; set subcell-indices
      (setf (getf properties 'SUBCELL-INDICES)
	    (loop+ (i (nr-dofs (subcell-ndofs vecfe)))
	       unless (zerop nr-dofs) collect i))
      ;; setup dofs
      (setq dofs
	    (loop+ (comp (fe components)
			 (local-off local-offset) (subcell-offset subcell-offsets))
	       nconcing
	       (loop+ ((dof (fe-dofs fe))) collecting
		  (new-vector-dof-from-dof dof comp subcell-offset))))
      )))

(defmethod make-local-vec ((vecfe <vector-fe>) &optional (multiplicity 1))
  (map 'simple-vector
       #'(lambda (fe) (make-real-matrix (nr-of-dofs fe) multiplicity))
       (components vecfe)))

(defmethod make-local-mat ((vecfe <vector-fe>))
  (let* ((n-comps (nr-of-components vecfe))
	 (result (make-array (list n-comps n-comps) :initial-element nil)))
    (loop+ (i (fe1 (components vecfe))) do
       (loop+ (j (fe2 (components vecfe))) do
	  (setf (aref result i j)
		(make-real-matrix (nr-of-dofs fe1) (nr-of-dofs fe2)))))
    result))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Interpolation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric interpolation-function (fe func &key &allow-other-keys)
  (:documentation "Assserts a correct form of @arg{func} for
interpolation with the finite element @arg{fe}.")
  (:method ((fe <scalar-fe>) func &key &allow-other-keys)
    "Asserts a scalar value of @arg{func}."
    #'(lambda (x)
	(let ((value (funcall func x)))
	(if (numberp value) value (aref value 0)))))
  (:method ((fe <vector-fe>) func &key dof)
    "Returns a function for the component of @arg{dof}."
    #'(lambda (x) (aref (funcall func x)
			(dof-component dof)))))

(defmethod interpolate-on-refcell ((fe <fe>) function)
  "Interpolates @arg{function} on the reference cell of the scalar finite
element @arg{fe}."
  (let ((values (loop for dof in (fe-dofs fe) when (interior-dof? dof)
		      collecting (evaluate dof (interpolation-function fe function :dof dof)))))
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
;;; <fe-discretization>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <fe-discretization> (<discretization>)
  ()
  (:documentation "FE discretization base class."))

(defgeneric get-fe (fe-disc cell)
  (:documentation "Returns the finite element for the given discretization
and reference cell."))

(defgeneric quadrature-rule (fe)
  (:documentation "Computes the quadrature rule to be used for the finite
element @arg{fe}."))

(defmethod quadrature-rule ((fe <fe>))
  "Standard quadrature rule for fe."
  (let ((refcell (reference-cell fe))
	(order (discretization-order fe)))
    (gauss-rule (mapcar #'dimension (factor-simplices refcell))
		;; does not integrate reaction terms precisely
		#+(or)(if (typep refcell '<simplex>) order (1+ order))
		(1+ order))))

(defclass <cell-fe-discretization> (<fe-discretization>)
  ((cell->fe :initarg :cell->fe :documentation
	     "A function mapping a cell to a finite element."))
  (:documentation "Finite element discretization where the finite elements
can differ from cell to cell.  Especially, hp-FEM are included."))

(defmethod get-fe ((disc <cell-fe-discretization>) cell)
  (funcall (slot-value disc 'cell->fe) cell))

(defclass <standard-fe-discretization> (<cell-fe-discretization>)
  ()
  (:documentation "Finite element discretization where the finite element
depends only on the type of the reference cell."))

(defmethod get-fe ((disc <standard-fe-discretization>) cell)
  (funcall (slot-value disc 'cell->fe)
	   (reference-cell cell)))

(defclass <scalar-fe-discretization> (<standard-fe-discretization>)
  ((order :reader discretization-order :initarg :order))
  (:documentation "Class for scalar fe discretizations."))

(defmethod nr-of-components ((fe-disc <scalar-fe-discretization>)) 1)

(defclass <vector-fe-discretization> (<standard-fe-discretization>)
  ((components :accessor components :initarg :components))
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
	(map 'vector
	     #'(lambda (comp-disc i)
		 (make-instance
		  '<scalar-fe-discretization>
		  :order (discretization-order comp-disc)
		  :cell->fe
		  #'(lambda (cell)
		      (aref (components (get-fe disc cell)) i))))
	     components (range< 0 (length components))))
  disc)


(defmethod discretization-order ((vecfe-disc <vector-fe-discretization>))
  (loop for disc across (components vecfe-disc)
	maximize (discretization-order disc)))

(defmethod nr-of-components ((vecfe-disc <vector-fe-discretization>))
  (length (components vecfe-disc)))

(defmethod component ((fedisc <fe-discretization>) i)
  (make-instance
   '<cell-fe-discretization> :cell->fe
   #'(lambda (cell)
       (aref (components (get-fe fedisc cell)) i))))

(defmethod component ((fedisc <vector-fe-discretization>) i)
  (aref (components fedisc) i))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ip-values, ip-gradients
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; The following two functions might be the standard way to generate the
;;; quadrature information given a finite element and a quadrature rule.
;;; They should be used in their memoized form

(defun base-function-values-at-ips (fe qrule)
  "Returns a list of nr-ip float-matrices of dimension (n-basis x 1)."
  (map 'vector
       #'(lambda (ip)
	   (make-real-matrix
	    (loop+ ((shape (fe-basis fe)))
	       collect (list (evaluate shape (ip-coords ip))))))
        (integration-points qrule)))
(memoize-symbol 'base-function-values-at-ips :test 'equal)

(defun base-function-gradients-at-ips (fe qrule)
  "Returns a list of nr-ip float-matrices of dimension (n-basis x dim)."
  (map 'vector
       #'(lambda (ip)
	   (make-real-matrix
	    (loop+ ((shape (fe-basis fe))) collect
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

(defmethod ip-values ((fe <scalar-fe>) qrule)
  "Returns a list of nr-ip float-matrices of dimension (n-basis x 1)."
  (funcall #'base-function-values-at-ips fe qrule))

(defmethod ip-gradients ((fe <scalar-fe>) qrule)
  (funcall #'base-function-gradients-at-ips fe qrule))

(defmethod ip-values ((fe <vector-fe>) qrule)
  "Returns a list of vectors of length components."
  (funcall #'vector-fe-ip-values fe qrule))

(defmethod ip-gradients ((fe <vector-fe>) qrule)
  (funcall #'vector-fe-ip-gradients fe qrule))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Coefficient input construction for fes
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun construct-coeff-input (cell global values gradients fe-parameters)
  "Constructs a coefficient input list from FE data @arg{cell} is the cell,
@arg{global} is the global coordinate of the integration point,
@arg{values} and @arg{gradients} the values and gradients of the shape
functions at the ip, and @arg{fe-parameters} are the corresponding data of
fe-functions to be evalutated."
  (list* :global global :cell cell
	 (loop for (obj data) on fe-parameters by #'cddr
	       collect (if (symbolp obj) obj (car obj))
	       collect
	       (let ((value (if (vectorp data)
				(map 'vector #'m*-tn data values)
				(m*-tn values data))))
		 (if (symbolp obj)
		     value
		     (list value
			   (if (vectorp data)
			       (map 'vector #'m*-tn data gradients)
			       (m*-tn gradients data))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Polynomial spaces on cells
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun P-nomials-of-degree (simplex deg &optional (type '=))
  "Returns a list of monomials of degree = deg or <= deg for a simplex."
  (n-variate-monomials-of-degree (dimension simplex) deg type))

(defun encapsulate (item dim)
  "Encapsulates @var{item} @var{dim} times."
  (if (zerop dim)
      item
      (encapsulate (list item) (- dim 1))))

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
  (let ((result vector)
	(coeffs (and basis (projection-coefficients vector basis pairing))))
    (loop+ ((coeff coeffs) (phi basis))
       doing (axpy! (- coeff) phi result))
    result))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; fe-cell-geometry
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun fe-cell-geometry (cell qrule &key metric volume)
  "Collects cell geometry information inside a property list."
  (loop for ip in (integration-points qrule)
	for local-coord = (ip-coords ip)
	for global-coord = (local->global cell local-coord)
	for Dphi = (local->Dglobal cell local-coord)
	for metric-ip = (and metric (funcall metric :local local-coord :global global-coord))
	for volume-ip = (and volume (funcall volume :local local-coord :global global-coord))
	for volume-at-point =
	(* (sqrt (abs (det (if metric
			       (m*-tn Dphi (m* metric-ip Dphi))
			       (m*-tn Dphi Dphi)))))
	   (or volume-ip 1.0))
	collect local-coord into local-coords
	collect global-coord into global-coords
	collect Dphi into gradients
	collect volume-at-point into volumes
	collect (when (= (nrows Dphi) (ncols Dphi)) (m/ Dphi)) into gradient-inverses
	collect (* (ip-weight ip) volume-at-point) into weights
	finally
	(return (list :cell cell :local-coords local-coords :global-coords global-coords
		      :gradients gradients :volume volumes
		      :gradient-inverses gradient-inverses :weights weights))))

;;; For the following to be effective, we should eliminate the consing by
;;; handing over a geometry to be filled
#+(or)
(defun fe-cell-geometry (cell qrule &key metric volume)
  "Collects cell geometry information inside a property list."
  (let* ((mapping (cell-mapping cell))
	 (fast-p (and (not metric) (not volume) (typep mapping '<linear-function>)))
	 Dphi origin volume-at-point Dphi-inverse)
    (when fast-p
      (setf origin (first (corners cell))
	    Dphi (evaluate-gradient mapping origin)
	    volume-at-point (* (sqrt (abs (det (m*-tn Dphi Dphi)))))
	    Dphi-inverse (when (= (nrows Dphi) (ncols Dphi)) (m/ Dphi))
	    ))
    (loop
     for ip in (integration-points qrule)
     for local-coord = (ip-coords ip)
     for local-weight = (ip-weight ip)
     for global-coord = (if fast-p
			    (gemm 1.0 Dphi local-coord 1.0 origin)
			    (evaluate mapping local-coord)) do
     (unless fast-p
       (break)
       (let ((metric-ip (and metric (funcall metric :local local-coord :global global-coord)))
	     (volume-ip (and volume (funcall volume :local local-coord :global global-coord))))
	 (setf Dphi (evaluate-gradient mapping local-coord)
	       Dphi-inverse (when (= (nrows Dphi) (ncols Dphi)) (m/ Dphi))	       
	       volume-at-point (* (sqrt (abs (det (if metric
						      (m*-tn Dphi (m* metric-ip Dphi))
						      (m*-tn Dphi Dphi)))))
				  (or volume-ip 1.0)))))
     collect local-coord into local-coords
     collect global-coord into global-coords
     collect Dphi into gradients
     collect volume-at-point into volumes
     collect Dphi-inverse into gradient-inverses
     collect (* local-weight volume-at-point) into weights
     finally
     (return (list :cell cell :local-coords local-coords :global-coords global-coords
		   :gradients gradients :volume volumes
		   :gradient-inverses gradient-inverses :weights weights)))))

;;; Testing:

(defun test-fe ()
  (Q-nomials-of-degree *unit-interval* 1 '=)
  (Q-nomials-of-degree *reference-vertex* 0 '<=)
  (Q-nomials-of-degree (n-cube 3) 2)
  #+(or) ; works only if lagrange.lisp has been evaluated
  (let ((gm (gram-matrix (Q-nomials-of-degree *unit-interval* 1 '<=)
			 (lagrange-dofs *unit-interval* 1) #'evaluate)))
    (assert (< (abs (- 1.0 (mref (m* (m/ gm) gm) 0 0)) ) 1.0e-10)))
  )

;;; (test-fe)
(fl.tests:adjoin-test 'test-fe)
