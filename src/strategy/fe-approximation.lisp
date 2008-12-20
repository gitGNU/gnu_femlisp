;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; fe-approximation.lisp - Solution strategy for FE approximation
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

;;;; This module provides a general pattern for iterative finite element
;;;; approximation which is used for interpolation with finite elements and
;;;; for the solution of stationary problems with finite elements.

(in-package :fl.strategy)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Data used
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defvar *eta-observe*
  (list "      ETA" "~9,2,2E"
	#'(lambda (blackboard)
	    (getbb blackboard :global-eta)))
  "Observe an estimate of the global error.")

(defvar *fe-approximation-observe*
  (list (list "Step" "~4D"
	      #'(lambda (blackboard) (getbb blackboard :step)))
	(list " CELLS" "~6D"
	      #'(lambda (blackboard)
		  (nr-of-surface-cells (getbb blackboard :mesh))))
	(list "    DOFS" "~8D"
	      #'(lambda (blackboard)
		  (with-items (&key solution) blackboard
		    (and solution (total-entries solution)))))
	*time-observe*)
  "Standard observe quantities for fe approximation.")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Class definition and interface
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <fe-approximation> (<strategy>)
  ((fl.iteration::observe :initform *fe-approximation-observe* :initarg :observe
    :documentation "Providing initform for <iteration> slot.")
   (plot-mesh :initform t :initarg :plot-mesh
    :documentation "Plot the mesh at the beginning and after changes.  Can
be a function in which case it is called on the mesh to do the plotting.")
   (fe-class :reader fe-class :initarg :fe-class
    :documentation "The class of finite element.  If it is not set, it is
automatically chosen.")
   (estimator :initform nil :initarg :estimator
    :documentation "The error estimator, which computes information on the error
distribution in a hash-table in the :ETA-field on the blackboard, as well
as a global estimate in :GLOBAL-ETA which can be used to terminate the
iteration.")
   (indicator :initform (make-instance '<uniform-refinement-indicator>)
	      :initarg :indicator :documentation
	      "The error indicator which marks cells for local refinement.
Usually, this procedure will be based on the error distribution
approximated in the :ETA-field on the blackboard."))
  (:documentation "This class describes iterative finite element
appoximation strategies."))

(defgeneric approximate (strategy blackboard)
  (:documentation "Compute a finite element approximation within the
setting given on the blackboard."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Iteration methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ensure-ansatz-space (blackboard)
  "Ensures a finite element ansatz space on the blackboard.  If a
discretization is not provided, it is selected depending on the problem.
If a mesh is not provided and if the domain is curved, we use precise cell
mappings if these are available.  Otherwise, an isoparametric approximation
is used with order @math{max(2,p+1)} where p is the discretization order.
This is reasonable, because the dual problem for estimators based on
duality is discretized with @math{p+1}, and because the refinement near the
boundary works significantly better with @math{p>=2}."
  (with-items (&key problem mesh ansatz-space fe-class base-level)
      blackboard
    (unless ansatz-space
      (assert problem)
      (ensure fe-class (select-discretization problem blackboard))
      (ensure mesh (let* ((domain (domain problem))
			  (parametric
			   (let ((chars (domain-characteristics domain)))
			     (and (getf chars :curved)
				  (if (getf chars :exact)
				      :from-domain
				      (lagrange-mapping
				       (max 2 (1+ (discretization-order fe-class)))))))))
		     (if (find-cell (rcurry #'typep '<boundary-cell>) domain)
			 (change-class (triangulate domain :parametric parametric) '<hierarchical-mesh>)
			 (uniformly-refined-hierarchical-mesh domain (or base-level 0)
							      :parametric parametric))))
      (ensure ansatz-space (make-fe-ansatz-space fe-class problem mesh)))
    ;; Ensure consistency of items on blackboard.  This might be dropped
    ;; later, because probably the ansatz-space should be considered the
    ;; fundamental approximation entity
    (ensure problem (problem ansatz-space))
    (ensure mesh (mesh ansatz-space))
    (ensure fe-class (fe-class ansatz-space))
    ))
    
(defmethod initially :after ((fe-strategy <fe-approximation>) blackboard)
  "Ensures a finite element ansatz space and a first approximation to the
solution on the blackboard."
  (with-items (&key ansatz-space solution fe-class)
      blackboard
    (when (slot-boundp fe-strategy 'fe-class)
      (setf fe-class (fe-class fe-strategy)))
    (ensure-ansatz-space blackboard)
    (ensure solution (choose-start-vector ansatz-space))
    ;; ensure data correctness and compatibility
    (check (mesh ansatz-space))
    ;;
    (with-slots (estimator indicator fe-class) fe-strategy
      ;; fe-class?
      (unless (slot-boundp fe-strategy 'estimator)
	(setq estimator (select-estimator (problem ansatz-space) blackboard)))
      (unless (slot-boundp fe-strategy 'indicator)
	(setq indicator (select-indicator (problem ansatz-space) blackboard))))
    ;;
    (whereas ((mesh-plotter (slot-value fe-strategy 'plot-mesh)))
      (if (functionp mesh-plotter)
	  (funcall mesh-plotter (mesh ansatz-space))
	  (plot (mesh ansatz-space))))))

(defmethod intermediate ((fe-strategy <fe-approximation>) blackboard)
  "Approximates and estimates the error for the current ansatz-space."
  (approximate fe-strategy blackboard)
  (whereas ((estimator (slot-value fe-strategy 'estimator)))
    (estimate estimator blackboard))
  (call-next-method))

(defmethod terminate-p ((fe-strategy <fe-approximation>) blackboard)
  "Sets some variables before passing control to the check of termination
criteria."
  (with-items (&key nr-levels max-level mesh)
      blackboard
    (setq nr-levels (nr-of-levels mesh))
    (setq max-level (top-level mesh))
  (call-next-method)))

(defmethod next-step ((fe-strategy <fe-approximation>) blackboard)
  "Enlarges the ansatz space."
  (with-items (&key mesh refinement-table refined-cells)
      blackboard
    (indicate (slot-value fe-strategy 'indicator) blackboard)
    (setf refined-cells
	  (let ((ht refinement-table))
	    (nth-value 1 (refine mesh :indicator #'(lambda (cell)
						     (nth-value 1 (gethash cell ht)))))))
    (when (get-property  (domain mesh) 'fl.mesh::extensible-p)
      (extend mesh :test #'(lambda (cell) (member-of-skeleton? cell refined-cells))))
    (whereas ((mesh-plotter (slot-value fe-strategy 'plot-mesh)))
      (if (functionp mesh-plotter) (funcall mesh-plotter mesh) (plot mesh)))
    (assert (not (skel-empty-p refined-cells))) ; otherwise the estimator should have approved
    ;; update solution and interpolation and projection operators
    (update-I-P-sol blackboard)))

