;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; strategy-cdr.lisp - Using the strategy interface
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Demos
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun laplace-1d-demo-a ()
  "Solve the equation -u''(x) = x^5, u(0)=u(1)=0 adaptively.  We use
fourth order finite elements, a direct solver, and error estimation by
projection between subsequent levels."
  (defparameter *result*
    (let* ((dim 1) (order 4)
	   (problem
	    (cdr-model-problem dim :source #'(lambda (x) #I"x[0]^^5"))))
      (solve
       (make-instance
	'<stationary-fe-strategy> :fe-class (lagrange-fe order)
	:estimator (make-instance '<projection-error-estimator>)
	:indicator (make-instance '<largest-eta-indicator> :fraction 0.5 :from-level 2)
	:success-if '(and (>= :nr-levels 2) (< :global-eta 1.0d-8))
	:solver *lu-solver*
	:output t)
       (blackboard :problem problem))))
  (sleep 1.0)
  (plot (getf *result* :solution)))

;;; (laplace-1d-demo-a)

(let ((demo
       (make-demo
	:name "1D-laplace-A" :short "u''=x^5, adaptive, projection error estimate"
	:long (documentation 'laplace-1d-demo-a 'function) 
	:execute #'adaptive-laplace-1d-a)))
  (adjoin-demo demo *laplace-demo*)
  (adjoin-demo demo *adaptivity-demo*))

(defun laplace-1d-demo-b ()
  "Solve the equation -u''(x) = x^5, u(0)=u(1)=0 adaptively.  We use third
order finite elements, a direct solver, and a duality error estimator for
the error in the load functional (which is the same as the energy error).
The dual problem is solved with a fifth order method so that this error
estimator is asymptotically exact.  Accuracy threshold is 1.0e-8."
  (defparameter *result*
    (let* ((dim 1) (order 4)
	   (problem (cdr-model-problem dim :source #'(lambda (x) #I"x[0]^^5"))))
      (solve
       (make-instance
	'<stationary-fe-strategy> :fe-class (lagrange-fe order)
	:estimator (make-instance '<duality-error-estimator> :functional :load-functional)
	:indicator (make-instance '<largest-eta-indicator> :from-level 1)
	:success-if '(< :global-eta 1.0d-10)
	:solver #+(or)(s1-reduction-amg-solver order :output t) #-(or)*lu-solver*)
       (blackboard :problem problem :output t)))))

;;; (laplace-1d-demo-b)
(let ((demo
       (make-demo
	:name "1D-laplace-B" :short "u''=x^5, adaptive, duality error estimate"
	:long (documentation 'laplace-1d-demo-b 'function)
	:execute #'laplace-1d-demo-b)))
  (adjoin-demo demo *laplace-demo*)
  (adjoin-demo demo *adaptivity-demo*))

(defun laplace-l-domain-demo (dim order threshold)
  "Solve the equation -Delta u(x) = 1 on the L-domain with Dirichlet b.c.
We use finite elements of order p, an AMG solver, and a duality error
estimator for the error in the load functional (which is the same as the
energy error).  The dual problem is solved with a (p+1)-method so that this
error estimator is asymptotically exact."
  (defparameter *result*
    (solve
     (make-instance
      '<stationary-fe-strategy> :fe-class (lagrange-fe order)
      :estimator
      (make-instance '<duality-error-estimator> :functional :load-functional)
      :indicator (make-instance '<largest-eta-indicator> :fraction 0.25)
      :success-if `(<= :global-eta ,threshold)
      :solver (s1-reduction-amg-solver order))
     (blackboard
      :problem (cdr-model-problem (l-domain dim)
				  :source (constant-coefficient 1.0))
      :output t)))
  (plot (getf *result* :solution)))
;;; (laplace-l-domain-demo)

(defun create-l-domain-demo (dim order threshold)
  (let ((demo
	 (make-demo
	  :name (format nil "~DD-L-domain" dim)
	  :short (format nil "Solution of Delta u = 1 on an L-domain in ~DD" dim)
	  :long (concatenate
		 'string
		 (documentation 'laplace-l-domain-demo 'function)
		 "Parameters: order=~D, threshold=~A")
	  :execute (lambda () (laplace-l-domain-demo dim order threshold)))))
    (adjoin-demo demo *laplace-demo*)
    (adjoin-demo demo *adaptivity-demo*)))

;;; activate the demos
(create-l-domain-demo 2 4 1.0e-5)
(create-l-domain-demo 3 1 1.0e-2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Testing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun strategy-cdr-tests ()
  
;;; duality error estimator

(laplace-1d-demo-a)

;; exakt: 1/12
;; Dualitaet 18, gesch 8.8e-4, wahr 2.95e-4
;; IP: gesch 9.18e-4, wahr 2.95e-4

(laplace-1d-demo-b)
(laplace-l-domain-demo)

;;; 1D test case: Dirichlet b.c.
(defparameter *result*
  (time
   (let* ((dim 1) (order 1)
	  (problem (cdr-model-problem
		    dim :source #'(lambda (x) #I"-4*exp(-2*x[0])")
		    :dirichlet
		    #'(lambda (x) #I"exp(-2*x[0])"))))
     (solve
      (make-instance
       '<stationary-fe-strategy> :fe-class (lagrange-fe order)
       :estimator (make-instance '<projection-error-estimator>)
       :indicator (make-instance '<largest-eta-indicator> :fraction 1.0)
       :success-if '(or (and (>= :nr-levels 2) (< :global-eta 1.0d-8))
		     (= :max-level 4))
       :solver #+(or)(s1-reduction-amg-solver order) #-(or)*lu-solver*)
      (blackboard :problem problem :output t)))))

)  ; end tests

(fl.tests:adjoin-test 'strategy-cdr-tests)