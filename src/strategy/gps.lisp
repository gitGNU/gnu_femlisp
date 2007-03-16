;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; gps.lisp - General problem solving strategy 
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

(in-package :fl.strategy)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Choice of a problem-dependent linear solver
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod select-linear-solver ((problem fl.cdr::<cdr-problem>) blackboard)
  (declare (ignore blackboard))
  (let ((dim (dimension (domain problem))))
    (make-instance
     '<safe-linear-solver>
     :iteration
     (if (and (= dim 1) (get-property problem 'coercive))
	 (make-instance '<lu> :store-p nil)
	 (let ((smoother (if (> dim 2) *gauss-seidel* (geometric-ssc))))
	   (make-instance '<s1-reduction> :max-depth 2
			  :smoother smoother :pre-steps 1 :post-steps 1
			  :gamma 2 :coarse-grid-iteration
			  (make-instance '<s1-coarse-grid-iterator>))))
     :success-if `(or (zerop :defnorm)
		   (and (>= :step 2) (> :step-reduction 0.9) (< :defnorm 1.0e-8)))
     :failure-if `(and (> :step 2) (> :step-reduction 0.9))
     )))

(defmethod select-linear-solver ((problem fl.elasticity::<elasticity-problem>) blackboard)
  (declare (ignore blackboard))
  (let ((dim (dimension (domain problem))))
    (make-instance
     '<safe-linear-solver> :iteration
     (let ((smoother (if (>= dim 3)
			 *gauss-seidel*
			 (geometric-ssc))))
       (geometric-cs :smoother smoother	:pre-steps 1 :post-steps 1 :gamma 2 :fmg t))
     :success-if `(or (zerop :defnorm)
		   (and (>= :step 2) (> :step-reduction 0.9) (< :defnorm 1.0e-8)
		     (> :reduction 1.0e-2)))
     :failure-if `(and (> :step 1) (> :step-reduction 0.9))
     )))

(defmethod select-linear-solver ((problem fl.navier-stokes::<navier-stokes-problem>) blackboard)
  (declare (ignore blackboard))
  (make-instance
   '<safe-linear-solver> :iteration
   (let ((smoother (make-instance '<vanka>)))
     (geometric-cs
      :coarse-grid-iteration
      ;;(make-instance '<lu>) #+(or)
      (make-instance '<multi-iteration> :nr-steps 1 :base smoother)
      :smoother smoother :pre-steps 1 :post-steps 1 :gamma 2))
     :success-if `(or (zerop :defnorm)
		   (and (>= :step 2) (> :step-reduction 0.9) (< :defnorm 1.0e-8)))
   :failure-if `(and (> :step 2) (> :step-reduction 0.9))
   ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Choice of problem-dependent error estimator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric select-estimator (problem blackboard)
  (:documentation "Select an error estimator."))

(defmethod select-estimator (problem blackboard)
  (declare (ignore problem))
  (with-items (&key functional) blackboard
    (when functional
      (make-instance '<duality-error-estimator> :functional functional))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Choice of problem-dependent error estimator
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric select-indicator (problem blackboard)
  (:documentation "Select a refinement indicator."))

(defmethod select-indicator ((problem t) blackboard)
  "For systems, there is still a problem with multigrid and local
refinements.  Therefore, we do not allow local refinement here."
  (declare (ignore blackboard))
  (make-instance '<uniform-refinement-indicator>))

(defmethod select-indicator ((problem <cdr-problem>) blackboard)
  "For CDR problems, solving with local refinement works.  If an error is
estimated, we refine those cells within some factor of the maximum error."
  (if (select-estimator problem blackboard)
      (make-instance '<largest-eta-indicator> :pivot-factor 0.1
		     :from-level 1 :block-p t)
      (make-instance '<uniform-refinement-indicator>)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; GPS = "General problem solver"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod solve ((blackboard blackboard) &optional dummy)
  "Implements the simplest interface for problem solving.  Only a
blackboard describing the problem is passed and the solution method is
open."
  (assert (null dummy))
  (with-items (&key strategy problem ansatz-space matrix rhs
		    fe-class solver estimator indicator)
      blackboard
    ;; the problem might be given already in a fixed ansatz-space
    (ensure problem (and ansatz-space (problem ansatz-space)))
    (cond
      (strategy				; there is a strategy on the blackboard
       (solve strategy blackboard))
      (problem				; there is a problem on the blackboard
       (setq
	strategy
	(typecase problem
	  (<lse> (select-linear-solver problem blackboard))
	  (<evp> (select-solver problem blackboard))
	  (<nlse> (select-solver problem blackboard))
	  (<pde-problem>
	   (apply #'make-instance '<stationary-fe-strategy>
		  (loop with flag = (cons nil nil)
			for item in '(:fe-class :solver :estimator :indicator
				      :success-if :failure-if :observe :plot-mesh
				      :base-level)
			for value = (getbb blackboard item flag)
			unless (eq value flag) collect item and collect value)))
	  (t (error "SOLVE does not know this problem."))))
       (let* ((bb-output (getbb blackboard :output :keep))
	      (*output-depth* (if (eq bb-output :keep)
				  *output-depth*
				  bb-output)))
	 (solve strategy blackboard)))
      ((and matrix rhs)			; there seems to be a linear problem on the blackboard
       (setq problem (lse :matrix matrix :rhs rhs))
       (solve blackboard))
      (t (error "SOLVE does not find a problem on the blackboard.")))))


