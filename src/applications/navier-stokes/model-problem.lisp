;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; model-problem.lisp - Model problems for Navier-Stokes
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

(in-package :application)

(defvar *stokes-demo*
  (make-demo :name "Stokes"
	     :short "Solve the Stokes equation on different domains"))
(adjoin-demo *stokes-demo* *equation-demo*)

(defun stokes-driven-cavity-demo (dim order levels &key plot output)
  "Solves the driven cavity problem for Stokes' equation."
  (defparameter *result*
    (solve
     (make-instance
      '<stationary-fe-strategy>
      :fe-class (navier-stokes-lagrange-fe order dim 1)
      :estimator (make-instance '<projection-error-estimator>)
      :indicator (make-instance '<largest-eta-indicator> :fraction 0.5 :from-level levels)
      :success-if `(= :nr-levels ,levels)
      :solver
      #+(or)
      (make-instance '<special-solver> :solver-function
		     #'(lambda (&key matrix rhs &allow-other-keys)
			 (m* (sparse-ldu matrix :ordering (boundary-last-order matrix))
			     rhs)))
      #-(or)  ; no local mg here yet
      (make-instance
       '<linear-solver> :iteration
       (let ((smoother (make-instance '<vanka>)))
	 (geometric-cs
	  :coarse-grid-iteration
	  (make-instance '<multi-iteration> :nr-steps 1 :base smoother)
	  :pre-steps 1 :pre-smooth smoother
	  :post-steps 1 :post-smooth smoother
	  :gamma 2))
       :success-if `(or (> :step 10) (< :defnorm 1.0e-10))
       :failure-if `(and (> :step 2) (> :step-reduction 0.9))
       :output (eq output :all))
      :output t)
     (blackboard :problem (driven-cavity dim) :base-level 0)))
  (when plot
    (let ((solution (getf *result* :solution)))
      ;; plot components of cell solution tensor
      (dotimes (i (1+ dim))
	(plot solution :component i :depth 2)
	(sleep 1.0)))))

;;; (stokes-driven-cavity-demo 2 2 1 :output :all :plot nil)
(let* ((dim 2) (order 3) (levels 4)
       (demo (make-demo
	      :name "driven-cavity"
	      :short "Solves the driven cavity problem for Stokes' equation"
	      :long
	      (format nil "~A~%~%Parameters: dim=~D, order=~D, levels=~D~%~%"
		      (documentation 'stokes-driven-cavity-demo 'function)
		      dim order levels)
	      :execute (lambda () (stokes-driven-cavity-demo dim order levels :plot t)))))
    (adjoin-demo demo *stokes-demo*))


;;;; Testing: (test-ns-model-problem)

(defun test-ns-model-problem ()
  
  (stokes-driven-cavity-demo 2 2 1 :output :all :plot nil)

  ;; Note: Solving the smooth DC-problem with order 1 FE does not work
  (defparameter *result*
    (time
     (let* ((dim 2) (level 2) (order 2) (delta 1)
	    (problem
	     (driven-cavity dim :smooth-p t)
	     #+(or)(periodic-cavity dim)
	     #+(or)
	     (standard-navier-stokes-problem
	      (?1 (n-cube-with-ellipsoidal-hole
		   2 :A (algebra::ellipse-matrix 0.25 0.3 0.7854))
		  (n-cube-domain dim))
	      :force (unit-vector-force dim 1))
	      )
	    (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	    (fedisc (navier-stokes-lagrange-fe order dim delta)))
       (multiple-value-bind (matrix rhs)
	   (discretize-globally problem h-mesh fedisc)
	 (let* ((smoother (make-instance '<vanka>))
		(cs (geometric-cs
		     :gamma 1 :base-level 1 :coarse-grid-iteration
		     (make-instance '<multi-iteration> :base smoother :nr-steps 4)
		     :pre-steps 1 :pre-smooth smoother
		     :post-steps 1 :post-smooth smoother))
		(ls (make-instance
		     '<linear-solver> :iteration cs
		     :success-if `(or (> :step 10) (< :defnorm 1.0e-10))
		     :failure-if `(and (> :step 2) (> :step-reduction 0.9))
		     :output t)))
	   (solve ls (blackboard :matrix matrix :rhs rhs))
	   ;;(linsolve matrix rhs :output t :iteration cs :maxsteps 10)
	   )))))
  (problem-info (driven-cavity 2 :smooth-p t))
  ;;(plot (getbb *result* :res) :component 0)
  ;;(plot (getbb *result* :solution) :component 2)
  ;;(plot (mesh (getbb *result* :res)))
  ;;(fe-extreme-values  (getbb *result* :solution))
)

(adjoin-femlisp-test 'test-ns-model-problem)
