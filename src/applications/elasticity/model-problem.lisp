;;; -*- mode: lisp; fill-column: 64 -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; model-problem.lisp - Model problems for linear elasticity
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

(defun elasticity-model-problem-computation (domain &key (output 1) plot)
  "Performs the model problem demo."
  (defparameter *result*
    (solve (blackboard
	    :problem (standard-elasticity-problem domain)
	    :plot-mesh t :output output :success-if `(> :time ,fl.demo:*demo-time*))))
  (when plot
    (plot (getbb *result* :solution))))
  
(defun make-elasticity-model-problem-demo (domain domain-name)
  (let ((title domain-name)
	(short (format nil "Solve an elasticity system on a ~A." domain-name))
	(long "Solves a linear elasticity problem with rhs
identical 1 on the given domain.  The solution strategy does
uniform refinement."))
    (let ((demo
	   (make-demo :name title :short short :long long
		      :execute (lambda ()
				 (elasticity-model-problem-computation domain :plot t)))))
      (adjoin-demo demo *elasticity-demo*))))

(make-elasticity-model-problem-demo (n-simplex-domain 2) "unit-triangle")
(make-elasticity-model-problem-demo (n-cube-domain 2) "unit-quadrangle")
(make-elasticity-model-problem-demo (n-simplex-domain 3) "unit-tetrahedron")
(make-elasticity-model-problem-demo (tensorial-domain '(1 2)) "unit-wedge-1-2")
(make-elasticity-model-problem-demo (tensorial-domain '(2 1)) "unit-wedge-2-1")
(make-elasticity-model-problem-demo (n-cube-domain 3) "unit-cube")

(defun test-elasticity-model-problem ()

;;; Linear elasticity problem
(defparameter *result*
  (time
   (let ((dim 2))
     (solve
      (blackboard
       :problem
       (standard-elasticity-problem
	(n-cube-domain dim) :lambda 1.0 :mu 1.0
	:force (constant-coefficient (make-array dim :initial-element (ones 1))))
       
;;        :fe-class (lagrange-fe 1 :nr-comps dim)
;;        :solver (make-instance '<linear-solver> :iteration (make-instance '<stueben>
;; 									 :cg-max-size 10 :output t)
;; 			      :success-if '(and (> :step 1) (< :defnorm 1.0e-10)))
       :output t :success-if '(> :nr-levels 2))))))
(plot (getbb *result* :solution) :component 0)

)

;;; (fl.application::test-elasticity-model-problem)
(fl.tests:adjoin-test 'test-elasticity-model-problem)
