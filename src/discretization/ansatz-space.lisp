;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ansatz-space.lisp - ansatz spaces
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2003-2006 Nicolas Neuss, University of Heidelberg.
;;; Copyright (C) 2007- Nicolas Neuss, University of Karlsruhe.
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
;;; <fe-discretization>
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <fe-discretization> (<discretization>)
  ()
  (:documentation "FE discretization base class."))

(defgeneric get-fe (fe-disc cell)
  (:documentation "Returns the finite element for the given discretization
and reference cell."))

(defclass <standard-fe-discretization> (<fe-discretization>)
  ((cell->fe :initarg :cell->fe :documentation
	     "A function mapping a cell to a finite element. Usually this
will be a closure memoized over the reference cell."))
  (:documentation "Finite element discretization where the finite element
depends only on the cell (usually via its reference cell)."))

(defmethod get-fe ((disc <standard-fe-discretization>) cell)
  (funcall (slot-value disc 'cell->fe) cell))

(defclass <scalar-fe-discretization> (<standard-fe-discretization>)
  ((order :reader discretization-order :initarg :order))
  (:documentation "Class for scalar fe discretizations."))

(defmethod nr-of-components ((fe-disc <scalar-fe-discretization>)) 1)

(defclass <vector-fe-discretization> (<standard-fe-discretization>)
  ((components :reader components :initarg :components))
  (:documentation "Vector FE discretization class."))

(defmethod initialize-instance :after ((disc <vector-fe-discretization>) &key &allow-other-keys)
  "If the slot @symbol{cell->fe} is unbound, it is derived from the
components vector."
  (unless (slot-boundp disc 'cell->fe)
    (setf (slot-value disc 'cell->fe)
	  (with-memoization (:id 'initialize-vector-fe-discretization :debug t)
	    (lambda (cell)
	      (memoizing-let ((refcell (reference-cell cell)))
                (break)  ; apparently not used at the moment
		(make-instance
		 '<vector-fe> :cell refcell :discretization disc :components
		 (map 'vector #'(lambda (comp-disc) (get-fe comp-disc refcell))
		      (components disc)))))))))

(defmethod discretization-order ((vecfe-disc <vector-fe-discretization>))
  (loop for disc across (components vecfe-disc)
	maximize (discretization-order disc)))

(defmethod nr-of-components ((vecfe-disc <vector-fe-discretization>))
  (length (components vecfe-disc)))

(defmethod component ((fedisc <vector-fe-discretization>) i)
  (aref (components fedisc) i))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ansatz-space
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space> (property-mixin)
  ((fe-class :reader fe-class :initarg :fe-class :type <fe-discretization> :documentation
	     "The finite element class for this ansatz space.")
   (problem :reader problem :initarg :problem :documentation
	    "The proplem for this ansatz space which determines essential constraints.")
   (mesh :reader mesh :initarg :mesh :type <mesh> :documentation
	    "The mesh for this ansatz space which determines hanging-node constraints."))
  (:documentation "A finite element ansatz space is determined by finite
element discretization, mesh and problem.  The constraints are stored in
the slot @var{properties}."))

(defmethod get-fe ((as <ansatz-space>) cell)
  (funcall (slot-value (fe-class as) 'cell->fe) cell))

(defgeneric hierarchical-mesh (aso)
  (:documentation "The hierarchical mesh for the given ansatz-space or ansatz-space object.")
  (:method ((as <ansatz-space>))
      (the <hierarchical-mesh> (values (mesh as)))))

(defmethod multiplicity ((as <ansatz-space>))
  (multiplicity (problem as)))

(defmethod nr-of-components ((as <ansatz-space>))
  (nr-of-components (fe-class as)))

(defmethod discretization-order ((as <ansatz-space>))
  (discretization-order (fe-class as)))

(defmethod initialize-instance :after ((as <ansatz-space>) &key &allow-other-keys)
  (call-hooks 'initialize-ansatz-space as))

(defun make-fe-ansatz-space (fe-class problem mesh)
  "Constructor of @class{<ansatz-space>}."
  (let ((n1 (nr-of-components fe-class))
	(n2 (nr-of-components problem)))
    (when (and n1 n2 (/= n1 n2))
      (error "Number of components of FE discretization and problem don't agree")))
  (make-instance '<ansatz-space>
                 :fe-class fe-class
                 :problem problem :mesh mesh))

(defgeneric set-constraints (ansatz-space)
  (:documentation "Computes the constraint matrices for this ansatz-space.
Constraints arise partially because of the discretization, e.g. hanging
nodes, and partially because of essential boundary conditions.  Of course,
these matrices change when mesh or discretization are adapted."))

(defun lagrange-ansatz-space (problem mesh &key (order 1) (type :uniform))
  "A constructor for a problem-dependent Lagrange fe.  Here, the number of
components may vary with the respective patch."
  (assert (eq (domain problem) (domain mesh)))
  (let ((disc (make-instance '<vector-fe-discretization>)))
    (setf (slot-value disc 'cell->fe)
	  (with-memoization ()
	    (lambda (cell)
	      (memoizing-let ((patch (patch-of-cell cell mesh))
			      (refcell (reference-cell cell)))
		(cell-lagrange-fe
		 refcell
		 (make-array (fl.problem::count-components (components-of-patch patch problem))
			     :initial-element order)
		 type)))))
    ;; return the generated ansatz-space
    (make-fe-ansatz-space disc problem mesh)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; ansatz-space objects
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <ansatz-space-object> (property-mixin)
  ((ansatz-space :reader ansatz-space :initarg :ansatz-space :type <ansatz-space>))
  (:documentation "Mixin for objects to which an ansatz-space is associated."))

(defmethod mesh ((aso <ansatz-space-object>))
  (mesh (ansatz-space aso)))

(defmethod hierarchical-mesh ((aso <ansatz-space-object>))
  (hierarchical-mesh (ansatz-space aso)))

(defmethod fe-class ((aso <ansatz-space-object>))
  (fe-class (ansatz-space aso)))

(defmethod problem ((aso <ansatz-space-object>))
  (problem (ansatz-space aso)))

(defmethod make-analog ((aso <ansatz-space-object>))
  (call-hooks 'make-analog-aso
              (make-instance (class-of aso) :ansatz-space (ansatz-space aso))
              aso))

