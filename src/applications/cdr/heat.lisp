;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; heat.lisp - Solve the time-dependent heat equation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2005 Nicolas Neuss, University of Heidelberg.
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

(defun heat-equation-computation (domain order levels &key (output 1) plot)
  "Performs the heat equation demo."
  (let* ((problem (cdr-model-problem
		   domain :initial (constant-coefficient 1.0)
		   :reaction (constant-coefficient 0.0)
		   :source (constant-coefficient #m(0.0))))
	 (rothe (make-instance
		 '<rothe> :model-time 0.0 :time-step 0.01
		 :stationary-success-if `(> :nr-levels ,levels)
		 :success-if '(>= :step 20)
		 :output output :plot plot)))
  (defparameter *result*
    (iterate rothe (blackboard
		    :problem problem :fe-class (lagrange-fe order)
		    :plot-mesh nil :output t)))))

(defun make-heat-equation-demo (domain domain-name)
  (let ((title domain-name)
	(short (format nil "Solve the heat equation on a ~A." domain-name))
	(long "Solve the heat equation with initial condition 1,
source 0 and zero boundary conditions on the given domain.
First, the initial value is projected L2-orthogonally on the
ansatz space.  Next, several steps of the Rothe method are
performed."))
    (adjoin-demo
     (make-demo
      :name title :short short :long long
      :execute
      (let* ((dim (dimension domain))
	     (order 4)
	     (levels (case dim (1 5) (2 2) (3 1))))
	(lambda ()
	  (heat-equation-computation domain order levels :plot t))))
      *heat-demo*)))

(make-heat-equation-demo (n-cube-domain 1) "unit-interval")
(make-heat-equation-demo (n-simplex-domain 2) "unit-triangle")
(make-heat-equation-demo (n-cube-domain 2) "unit-quadrangle")

(defun test-heat-equation ()
  (heat-equation-computation (n-cube-domain 1) 4 5 :plot t)
  (heat-equation-computation (n-cube-domain 2) 3 2 :plot t)
  )

;;; (fl.application::test-heat-equation)
(fl.tests:adjoin-test 'test-heat-equation)
