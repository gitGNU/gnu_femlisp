;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; mg-cdr.lisp - Solving CDR problems with multigrid
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

(time
 (let* ((dim 3) (level 2) (order 4)
	(problem (cdr-model-problem dim))
	(smoother #+(or)(geometric-ssc)
		  #-(or)(make-instance '<gauss-seidel>))
	(cs (geometric-cs
	     :gamma 1 :pre-steps 1 :pre-smooth smoother :post-steps 0
	     :base-level (1- level)))
	(solver (make-instance '<linear-solver> :iteration cs
			       :success-if '(or (< :defnorm 1.0e-12) (> :step 20))))
	(mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	(fedisc (lagrange-fe order)))
   (multiple-value-bind (A b)
       (discretize-globally problem mesh fedisc)
     (solve solver (blackboard :matrix A :rhs b)))))

(time
 (let* ((dim 2) (level 4)
	(problem
	 (cdr-model-problem
	  dim :source #'(lambda (x) (if (>= (aref x 0) 0.5) 1.0d0 -1.0d0))))
	(v-cycle (geometric-cs :fmg t :base-level 1 :coarse-grid-iteration
			       (make-instance '<multi-iteration>
					      :base *gauss-seidel* :nr-steps 1))))
   (multiple-value-bind (A b)
      (problem-discretization problem :level level :order 1)
    (let ((sol (linsolve A b :output t :iteration v-cycle :maxsteps 2)))
      (plot sol)
      ))))

;;; geometric V-cycle for higher-order problems
(let* ((dim 2) (level 3) (order 1)
       (problem (cdr-model-problem dim))
       (v-cycle (geometric-cs :base-level 1)))
  (multiple-value-bind (A b)
      (problem-discretization problem :level level :order order)
    #+(or)
    (let ((mg-data (multilevel-decomposition v-cycle A)))
      (show (aref (getbb mg-data :a-vec) 0)))
    #+(or)
    (let ((mg-data (multilevel-decomposition v-cycle A)))
      (show (aref (getbb mg-data :i-vec) 0)))
    #+(or)(plot (mesh b))
    #-(or)
    (setq *result*
	  (linsolve A b :output t :iteration v-cycle
		    :maxsteps 10 :threshold 1.0e-10))
    ))
(plot *result*)

(defun mg-cdr-tests ()
;;; test if CR is small enough
(let* ((dim 2) (level 3) (order 1)
       (problem (cdr-model-problem dim)))
  (multiple-value-bind (A b)
      (problem-discretization problem :level level :order order)
    (let ((geomg (geometric-cs :base-level 1)))
      (let ((solver (make-instance '<linear-solver> :iteration geomg
				   :success-if '(> :step 10))))
	(setq *result* (solve solver (blackboard :matrix A :rhs b)))
	(let ((cr (getf (getbb *result* :report) :convergence-rate)))
	  (assert (< cr 0.08)))
	  (plot (getbb *result* :solution))
	  *result*))))
)
