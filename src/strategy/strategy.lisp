;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; strategy.lisp - Solution strategies
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

;;;; An optimal solution requires the combination of discretization, error
;;;; estimators, grid refinement and fast solution algorithms.  This file
;;;; provides some infrastructure for this combination.

(defgeneric solve-with (strategy problem &rest rest &key output &allow-other-keys)
  (:documentation "This generic function depends directly on strategy and
problem.  But usually, solve will be used instead.  The list of parameters
should be modifiable."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Standard interface for iterative problem solving
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <strategy> ()
  ((output :reader output :initarg :output :initform nil))
  (:documentation "General abstract problem solver class.  The interface of
problem solving is given by generic functions which act on an assembly
line."))

(defgeneric initial-guess (strategy assembly-line)
  (:documentation "Generates an initial guess for a solution to the
problem.  Pre-condition: problem.  Post-condition: solution."))

(defgeneric sufficient-p (strategy assembly-line)
  (:documentation "Works on the assembly line and sets end-p, if some end
criterion is met."))

(defgeneric improve-guess (strategy assembly-line)
  (:documentation "Improves the approximation to the solution.
Pre-condition: solution, defect.  Post-condition: updated solution and
defect."))

(defparameter *standard-strategy-observe-quantities*
  (list
   (list :global-eta " ETA~8T"
	 #'(lambda (assembly-line)
	     (with-items (&key global-eta) assembly-line
	       (if global-eta
		   (format t "~12,5,2E" global-eta)
		   (format t "~12T")))))
   (list :load-functional
	 "   LOAD FUNCTIONAL"
	 #'(lambda (assembly-line)
	     (with-items (&key solution rhs) assembly-line
	       (if (and solution rhs)
		   (format t "  ~19,12,2E" (dot solution rhs))
		   (format t "~21T"))))))
  "This is a list of possible entries for the :observe entry in
*strategy-output*.")

(defparameter *strategy-output*
  (let ((start-time nil))
    (list
     :initial
     #'(lambda (assembly-line)
	 (declare (ignore assembly-line))
	 (format t "~& CELLS     DOFS  MENTRIES   TIME") ; 35 characters
	 (dolist (quantity (getf *strategy-output* :observe))
	   (format t (second quantity)))
	 (terpri)
	 (setq start-time (get-internal-run-time)))
     :after-step
     #'(lambda (assembly-line)
	 (with-items (&key mesh matrix solution accurate-p)
	   assembly-line
	   (when accurate-p
	     (format t "~&~5D~10D~10D~7,1F"
		     (nr-of-surface-cells mesh)
		     (and solution (total-entries solution))
		     (and matrix (total-entries matrix))
		     (float (/ (- (get-internal-run-time) start-time)
			       internal-time-units-per-second)))
	     (dolist (quantity (getf *strategy-output* :observe))
	       (funcall (third quantity) assembly-line))
	     (terpri))))
     :plot-mesh t :observe ())))

(defun strategy-observe (quantity)
  "Adds quantity to the observed quantities."
  (aif (assoc (car quantity) (getf *strategy-output* :observe))
       (setf (cdr it) (cdr quantity))
       (push quantity (getf *strategy-output* :observe))))

;;; observe the error estimator
(strategy-observe (assoc :global-eta *standard-strategy-observe-quantities*))

(defmethod solve-with ((strategy <strategy>) (problem <problem>)
		       &rest parameters &key (output *strategy-output*))
  (let ((assembly-line (apply #'make-assembly-line parameters)))
    (setf (get-al assembly-line :strategy) strategy)
    (setf (get-al assembly-line :problem) problem)
    (aand output (getf *strategy-output* :initial)
	  (funcall it assembly-line))
    (initial-guess strategy assembly-line)
    (loop (sufficient-p strategy assembly-line)
	  (aand output (getf *strategy-output* :after-step)
		(funcall it assembly-line))
	  (awhen (get-al assembly-line :after-step)
	    (funcall it assembly-line))
	  (when (get-al assembly-line :end-p)
	    (aand output (getf *strategy-output* :final)
		  (funcall it assembly-line))
	    (return assembly-line))
	  (improve-guess strategy assembly-line))))

(defmethod solve ((strategy <strategy>) &rest parameters &key problem &allow-other-keys)
  "The usual interface for problem solving."
  (apply #'solve-with strategy problem parameters))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Finite element problem solving
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <fe-strategy> (<strategy>)
  ((fe-class :reader fe-class :initform (ext:required-argument) :initarg :fe-class)
   (estimator :reader estimator :initform nil :initarg :estimator)
   (indicator :reader indicator :initform nil :initarg :indicator)
   (appraise :reader appraise :initarg :appraise :type function)
   (solver :reader solver :initarg :solver :type <solver>))
  (:documentation "Solution of problems by finite element methods. Here,
initial-guess puts a mesh, an ansatz-space and a guess for the solution on
the assembly line.  Also improve-guess is further specified to do an
error-estimation, refinement indication, refinement, new assembly, and
computation of a new solution.

Ideally, also the fe-class should be determined automatically inside the
cycle.  But we are not there, yet."))


;;; initial-guess

(defmethod initial-guess ((fe-strategy <fe-strategy>) assembly-line)
  "Constructs a finite element ansatz space and a first approximation to
the solution.  If a mesh is not provided on the assembly line and if the
domain is curved, we use precise cell mappings if these are available.
Otherwise, an isoparametric approximation is used with order p+1.  This is
reasonable, first, because the dual problem for the duality-error-estimator
is discretized with p+1, and; second, because the refinement of boundary
cells works better with p>=2."
  (with-items (&key problem mesh ansatz-space solution
		    assembled-p accurate-p base-level)
    assembly-line
    (let ((fe-class (fe-class fe-strategy))
	  (domain (domain problem)))
      ;; 
      (setf mesh (or mesh
		     (uniformly-refined-hierarchical-mesh
		      domain (or base-level 0) :parametric
		      (let ((chars (domain-characteristics domain)))
			(and (getf chars :curved)
			     (if (getf chars :exact)
				 :from-domain
				 (lagrange-mapping (1+ (discretization-order fe-class)))))))))
      (awhen (getf *strategy-output* :plot-mesh)
	(if (functionp it) (funcall it mesh) (plot mesh)))
      (setf ansatz-space (make-fe-ansatz-space fe-class problem mesh))
      (setf solution (make-ansatz-space-vector ansatz-space))
      (setf assembled-p nil accurate-p nil))))

;;; sufficient-p

(defmethod sufficient-p ((fe-strategy <fe-strategy>) assembly-line)
  (with-items (&key assembled-p accurate-p end-p)
    assembly-line
    (when (and assembled-p accurate-p)
      (awhen (estimator fe-strategy)
	(estimate it assembly-line))
      (setf end-p (apply (appraise fe-strategy) assembly-line)))))

;;; improve-guess

(defmethod improve-guess ((fe-strategy <fe-strategy>) assembly-line)
  (with-items (&key mesh ansatz-space assembled-p accurate-p
		    interpolation projection solution interior-matrix
		    interior-rhs refinement-table refined-cells)
      assembly-line
    (with-slots (solver indicator)
      fe-strategy
      ;;
      (when accurate-p
	;; the solution on this ansatz-space was accurate: we improve the
	;; approximation by enlarging the ansatz-space
	(indicate indicator assembly-line)
	(setf refined-cells
	      (let ((ht refinement-table))
		(nth-value 1 (refine mesh :test #'(lambda (cell) (gethash cell ht))))))
	(when (extensible-p (domain mesh))
	  (extend mesh :test #'(lambda (cell) (member-of-skeleton? cell refined-cells))))
	(awhen (getf *strategy-output* :plot-mesh)
	  (if (functionp it) (funcall it mesh) (plot mesh)))
	(assert (not (skel-empty-p refined-cells))) ; otherwise the estimator should have approved
	;; update solution and interpolation and projection operators
	(update-I-P-sol assembly-line)
	(setf assembled-p nil accurate-p nil))
      ;;
      (unless assembled-p
	;; assemble if necessary (better would be a local assembly)
	(setf interior-matrix nil interior-rhs nil)
	(fe-discretize assembly-line)
	(setf assembled-p t)
	(setf accurate-p nil))
      ;;
      (unless accurate-p
	;; improve approximation by solving better
	(setf solution (apply #'solve solver assembly-line))
	;; We assume for the moment that one step was sufficient (direct
	;; solver or good iterative scheme).  One could also leave
	;; accurate-p nil.  In this case, an additional convergence test
	;; should check for convergence in this ansatz-space.
	(setf accurate-p t))
      )))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Standard routine for appraise
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun stop-if (&key nr-levels>= eta<=)
  #'(lambda (&key assembled-p accurate-p mesh global-eta &allow-other-keys)
      (and assembled-p accurate-p
	   (or (and nr-levels>= (>= (nr-of-levels mesh) nr-levels>=))
	       (and eta<= global-eta (<= global-eta eta<=))))))

;;; appraise

(defun global-estimate-smaller-than (threshold)
  "Checks if global-eta is smaller than the given threshold."
  #'(lambda (&key assembled-p accurate-p global-eta &allow-other-keys)
      (and assembled-p accurate-p global-eta (< global-eta threshold))))

(defun nr-of-levels>= (level)
  "Checks if the given number of levels is reached."
  #'(lambda (&key mesh assembled-p accurate-p &allow-other-keys)
      (and assembled-p accurate-p (>= (nr-of-levels mesh) level))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; some strategies
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun linear-problem-fe-strategy (&key (order 1) threshold)
  (make-instance
   '<fe-strategy> :fe-class (lagrange-fe order)
   :estimator (make-instance '<projection-error-estimator>)
   :indicator (make-instance '<largest-eta-indicator>)
   :appraise (global-estimate-smaller-than threshold)
   :solver *lu-solver*))

