;;; -*- mode: lisp; fill-column: 64 -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; amg-cdr.lisp - Solving CDR problems with AMG
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

(in-package :fl.application)

(defvar *amg-cdr-demo*
  (make-demo :name "amg-cdr"
	     :short "AMG solving of C-D-R eqns"))
(adjoin-demo *amg-cdr-demo* *laplace-demo*)

(defun amg-cdr-solver ()
  "The standard AMG solver for the demos in this file."
  (make-instance '<linear-solver>
		 :iteration (make-instance '<stueben> :cg-max-size 16 :output t)
		 :success-if '(or (< :defnorm 1.0e-10)
			       (and (> :step 2) (> :step-reduction 0.9)))
		 :failure-if '(> :time 10.0)))

;;; A demo 

(adjoin-demo
 (make-demo
  :name "anisotropy-Q" :short "Solving anisotropy with Q1-FE"
  :long
  "Solves $-eps u_xx - u_yy = 1$ discretized with Q1-FE on
aligned grids."
  :execute
  (lambda ()
    (let ((eps (float (user-input "eps: " #'realp) 1.0)))
      (defparameter *result*
	(solve (blackboard
		:problem
		(cdr-model-problem
		 2 :diffusion (constant-coefficient (diag (double-vec eps 1.0))))
		:output t :plot-mesh t
		:success-if `(> :time 10.0)
		:fe-class (lagrange-fe 1) :base-level 1
		:solver (amg-cdr-solver))))))
  :test-input "0.01~%")
 *amg-cdr-demo*)

(defun test-amg-cdr ()

;;; Testing the S1-reduction AMG
(time
 (let* ((dim 2) (level 4) (order 2)
	(problem (cdr-model-problem dim)))
   (multiple-value-bind (A b)
       (problem-discretization problem :level level :order order)
     (let ((amg (s1-reduction-amg-solver order :reduction 1.0d-10 :output t)))
       (solve (blackboard :solver amg :matrix A :rhs b :output t))))))

;;; k-fold jump in refinement depth
(time
 (let* ((dim 2) (order 2)
	(problem (cdr-model-problem dim))
	(fedisc (lagrange-fe order))
	(h-mesh (make-hierarchical-mesh-from-domain (domain problem))))
   #-(or) (loop repeat 1 do (refine h-mesh))
   #-(or) (loop repeat 2 do
	       (refine h-mesh :indicator (rcurry #'inside-cell? (make-double-vec dim 0.25))))
   #-(or) (plot h-mesh)
   #-(or)
   (multiple-value-bind (mat rhs)
       (discretize-globally problem h-mesh fedisc)
     (let ((amg (s1-reduction-amg-solver order :maxsteps 10 :reduction 1.0e-10)))
       (solve (blackboard :solver amg :matrix mat :rhs rhs :output t))))))

(time
 (let* ((dim 2) (order 2)
	(problem (cdr-model-problem dim))
	(fedisc (lagrange-fe order))
	(h-mesh (make-hierarchical-mesh-from-domain (domain problem))))
   #-(or) (loop repeat 1 do (refine h-mesh))
   #-(or) (loop repeat 2 do
	       (refine h-mesh :indicator (rcurry #'inside-cell? #d(0.25 0.25))))
   #+(or) (plot h-mesh)
   #-(or)
   (multiple-value-bind (mat rhs)
       (discretize-globally problem h-mesh fedisc)
     (let ((amg (s1-reduction-amg-solver order :reduction 1.0e-5 :maxsteps 5)))
       (solve (blackboard :solver amg :matrix mat :rhs rhs :output t))))))
)

;;; (test-amg-cdr)
(fl.tests:adjoin-test 'test-amg-cdr)