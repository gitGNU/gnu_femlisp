;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; driven-cavity.lisp - Driven cavity computations
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

(defun ns-driven-cavity-demo (dim order levels &key plot output (reynolds 0.0))
  (defparameter *result*
    (solve 
     (blackboard
      :fe-class (navier-stokes-lagrange-fe order dim 1)
      :problem (driven-cavity dim :reynolds reynolds) :base-level 0
      :success-if (if levels
		      `(= :nr-levels ,levels)
		      `(> :time ,*demo-time*))
      :output output :observe
      (append *stationary-fe-strategy-observe*
	      (list
	       (list (format nil "~{                 u~1D~}" (range 1 dim))
		     "~{~19,10,2E~}"
		     #'(lambda (blackboard)
			 (let ((val (fe-value (getbb blackboard :solution)
					      (make-double-vec dim 0.5))))
			   (loop for i below dim collect (vref (aref val i) 0))))))))))
  (when plot
    (let ((solution (getbb *result* :solution)))
      ;; plot components of cell solution tensor
      (dotimes (i (1+ dim))
	(plot solution :component i :depth 2)
	(sleep 1.0)))))

;;; (ns-driven-cavity-demo 2 2 3 :output :all :plot nil :reynolds 100.0)
;;; (plot (getbb *result* :solution) :component 0 :depth 2)

(defun make-driven-cavity-demo (dim order reynolds)
  "DC-~dim~D-~reynolds~ - Solves the driven cavity problem (Re=~reynolds~).

Solve the ~dim~D driven cavity problem for the Navier-Stokes equation using
Taylor-Hood finite elements (Q^{k+1})^~dim~/Q^k with k=~order~."
  (multiple-value-bind (title short long)
      (extract-demo-strings
       (documentation 'make-driven-cavity-demo 'function)
       `(("~dim~" . ,dim) ("~order~" . ,order) ("~reynolds~" . ,reynolds)))
    (let ((demo
	   (make-demo
	    :name title :short short :long long :execute
	    (lambda ()
	      (ns-driven-cavity-demo dim order 3 :output 1 :plot t
				     :reynolds (float reynolds 1.0))))))
      (adjoin-demo demo *navier-stokes-demo*))))

;;(ns-driven-cavity-demo 2 2 3 :output :all :plot t :reynolds 100.0)
(make-driven-cavity-demo 2 2 0)
(make-driven-cavity-demo 2 2 100)

;;;; Testing:

(defun test-driven-cavity ()
  
  (ns-driven-cavity-demo 2 2 1 :output :all :plot nil)
  (let ((sol (getbb *result* :solution)))
    (fe-value sol #d(0.5 0.5)))
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
	   (solve ls (blackboard :problem (lse :matrix matrix :rhs rhs)))
	   ;;(linsolve matrix rhs :output t :iteration cs :maxsteps 10)
	   )))))
  (describe (driven-cavity 2))
  ;;(plot (getbb *result* :res) :component 0)
  ;;(plot (getbb *result* :solution) :component 2)
  ;;(plot (mesh (getbb *result* :res)))
  ;;(fe-extreme-values  (getbb *result* :solution))
)

;;; (application::test-driven-cavity)
(fl.tests:adjoin-test 'test-driven-cavity)
