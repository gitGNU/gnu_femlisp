;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; stationary.lisp - Solution strategy for stationary problems
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

(in-package :strategy)

(defclass <strategy> (<iteration>)
  ()
  (:documentation "Strategies are iterations devoted to solving continuous
problems."))
  
(defgeneric solve-with (strategy problem blackboard)
  (:documentation "This generic function can be dispatched on both strategy
and problem."))

(defmethod solve-with ((strategy <strategy>) (problem <problem>) blackboard)
  (setf (getbb blackboard :problem) problem)
  (solve strategy blackboard))

(defmethod solve ((strategy <strategy>) blackboard)
  "The usual interface for problem solving."
  (setf (getbb blackboard :strategy) strategy)
  (iterate strategy blackboard))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Finite element solving of stationary problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *stationary-fe-strategy-observe*
  (list (list "Step" "~4D"
	      #'(lambda (blackboard) (getbb blackboard :step)))
	(list " CELLS" "~6D"
	      #'(lambda (blackboard)
		  (nr-of-surface-cells (getbb blackboard :mesh))))
	(list "      DOFS" "~10D"
	      #'(lambda (blackboard)
		  (with-items (&key matrix solution) blackboard
		    (and matrix solution
			 (* (total-nrows matrix) (multiplicity solution))))))
	(list "  MENTRIES" "~10D"
	 #'(lambda (blackboard)
	     (with-items (&key matrix) blackboard
	       (and matrix (total-entries matrix)))))
	;; inside a closure
	(let ((start-time))
	  (list #'(lambda ()
		    (setq start-time (get-internal-run-time))
		    "   TIME")
		"~7,1F"
		#'(lambda (blackboard)
		    (declare (ignore blackboard))
		    (float (/ (- (get-internal-run-time) start-time)
			      internal-time-units-per-second))))))
  "Standard observe quantities for stationary fe-strategies.")

(defparameter *eta-observe*
  (list "         ETA" "~12,2,2E"
	#'(lambda (blackboard)
	    (getbb blackboard :global-eta))))

(defclass <stationary-fe-strategy> (<strategy>)
  ((iteration::observe
    :initform *stationary-fe-strategy-observe*
    :documentation "Providing initform for <iteration> slot.")
   (plot-mesh :initform t :initarg :plot-mesh
	      :documentation "Plot mesh at the beginning and after changes.
Can be a function in which case it is called on the mesh to do the
plotting.")
   (fe-class :reader fe-class :initform (ext:required-argument) :initarg :fe-class
	     :documentation "The class of fe.  Later on, this should be
automatically determined inside an hp-method.")
   (estimator :initform nil :initarg :estimator :documentation
	      "The error estimator, which computes information on the error
distribution in a hash-table in the :ETA-field on the blackboard, as well
as a global estimate in :GLOBAL-ETA which can be used to terminate the
iteration.")
   (indicator :initform nil :initarg :indicator :documentation
	      "The error indicator which marks cells for local refinement.
Usually, this procedure will be based on the error distribution
approximated in the :ETA-field on the blackboard.")
   (solver :initarg :solver :documentation
	   "The solver to be used for solving the discretized systems."))
  (:documentation "This class describes some iterative finite element
solution strategies for continuous, stationary PDE problems."))

(defmethod initially ((fe-strategy <stationary-fe-strategy>) blackboard)
  "Constructs a finite element ansatz space and a first approximation to
the solution.  If a mesh is not provided on the blackboard and if the
domain is curved, we use precise cell mappings if these are available.
Otherwise, an isoparametric approximation is used with order p+1.  This is
reasonable, first, because the dual problem for the duality-error-estimator
is discretized with p+1, and; second, because the refinement near the
boundary works significantly better with p>=2."
  (with-items (&key problem mesh ansatz-space solution base-level)
    blackboard
    (let ((fe-class (fe-class fe-strategy))
	  (domain (domain problem)))
      (setf mesh (or mesh
		     (uniformly-refined-hierarchical-mesh
		      domain (or base-level 0) :parametric
		      (let ((chars (domain-characteristics domain)))
			(and (getf chars :curved)
			     (if (getf chars :exact)
				 :from-domain
				 (lagrange-mapping (1+ (discretization-order fe-class)))))))))
      (check mesh)
      (whereas ((mesh-plotter (slot-value fe-strategy 'plot-mesh)))
	(if (functionp mesh-plotter) (funcall mesh-plotter mesh) (plot mesh)))
      (setf ansatz-space (make-fe-ansatz-space fe-class problem mesh))
      (setf solution (make-ansatz-space-vector ansatz-space))))
  (call-next-method))

(defmethod intermediate ((fe-strategy <stationary-fe-strategy>) blackboard)
  "Ensures accuracy of the solution and the error estimate."
  (with-items (&key mesh ansatz-space solver-blackboard 
		    interpolation projection solution rhs matrix
		    interior-matrix interior-rhs refinement-table refined-cells)
      blackboard
    (with-slots (solver estimator indicator) fe-strategy
      ;; assemble (better would be a local assembly)
      (setf interior-matrix nil interior-rhs nil)
      (fe-discretize blackboard)
      ;; improve approximation by solving
      (setf solver-blackboard
	    (solve solver (blackboard :matrix matrix :rhs rhs :solution solution)))
      (setf solution (getbb solver-blackboard :solution))
      ;; estimate
      (whereas ((estimator (slot-value fe-strategy 'estimator)))
	(estimate estimator blackboard))))
  (call-next-method))
  
(defmethod terminate-p ((fe-strategy <stationary-fe-strategy>) blackboard)
  "When there is an accurate solution we call the error estimator, if
provided, before passing control to the check of termination criteria."
  (with-items (&key nr-levels max-level mesh)
      blackboard
    (setq nr-levels (nr-of-levels mesh))
    (setq max-level (top-level mesh))
  (call-next-method)))

(defmethod next-step ((fe-strategy <stationary-fe-strategy>) blackboard)
  "Enlarges the ansatz space."
  (with-items (&key mesh refinement-table refined-cells)
      blackboard
    (indicate (slot-value fe-strategy 'indicator) blackboard)
    (setf refined-cells
	  (let ((ht refinement-table))
	    (nth-value 1 (refine mesh :test #'(lambda (cell) (gethash cell ht))))))
    (when (extensible-p (domain mesh))
      (extend mesh :test #'(lambda (cell) (member-of-skeleton? cell refined-cells))))
    (whereas ((mesh-plotter (slot-value fe-strategy 'plot-mesh)))
      (if (functionp mesh-plotter) (funcall mesh-plotter mesh) (plot mesh)))
    (assert (not (skel-empty-p refined-cells))) ; otherwise the estimator should have approved
    ;; update solution and interpolation and projection operators
    (update-I-P-sol blackboard)))

