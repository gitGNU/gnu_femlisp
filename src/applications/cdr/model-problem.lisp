;;; -*- mode: lisp; -*-

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

;;; exact values for u(1/2,...,1/2) for the solution of
;;;    - \delta u = 1
;;; in the unit square with Dirichlet boundary conditions

(defun exact-1d-solution ()
  "The precise value for x=1/2 is easily seen to be 1/8.
On the other hand it is equal to
1/pi^2 (4/pi)^d sum_{k odd} +/- 1/k^3"
  (let ((N 1000))
    (* (/ 4.0 (expt pi 3))
       (loop for i1 from N downto 0
	     for sign = (if (evenp i1) 1 -1) then (- sign)
	     for k1 = (float (1+ (* 2 i1)) 1.0)
	     sum (/ sign (* k1 k1 k1))))))

(defun exact-2d-solution ()
  "With the help of Fourier series we obtain the representation
1/pi^2 (4/pi)^d sum_{k_i odd} +/- 1/{k_1 k_2 (k_1^2+k_2^2)}
for the value u(1/2,1/2).
+/- means sin(k_1*pi/2) * sin (k_2*pi/2)"
  (time
   (let* ((N 1000)
	  (sign0 (expt -1.0 (if (evenp N) 0 1))))
     (* (/ 16.0 (expt pi 4))
	(loop for i1 from N downto 0
	      for sign1 double-float = sign0 then (- sign1)
	      for k1 double-float = (float (1+ (* 2 i1)) 1.0) summing
	      (loop for i2 from N downto 0
		    for sign = (* sign0 sign1) then (- sign)
		    for k2 double-float = (float (1+ (* 2 i2)) 1.0)
		    sum (/ sign
			   (* k1 k2 (+ (* k1 k1) (* k2 k2))))))))))

(defun exact-3d-solution ()
  "With the help of Fourier series we obtain the representation
1/pi^2 (4/pi)^d sum_{k_i odd} +/- 1/{k_1 k_2 k_3 (k_1^2+k_2^2+k_3^2)}
for the value u(1/2,1/2,1/2).
+/- means sin(k_1*pi/2) * sin (k_2*pi/2) * sin (k_3*pi/2)"
  (time
   (let* ((N 102)
	  (sign0 (expt -1.0 (if (evenp N) 0 1))))
     (* (/ 64.0 (expt pi 5))
	(loop for i1 from N downto 0
	      for sign1 = sign0 then (- sign1)
	      for k1 = (float (1+ (* 2 i1)) 1.0) summing
	      (loop for i2 from N downto 0
		    for sign2 = (* sign0 sign1) then (- sign2)
		    for k2 = (float (1+ (* 2 i2)) 1.0) summing
		    (loop for i3 from N downto 0
			  for sign = (* sign0 sign2) then (- sign)
			  for k3 = (float (1+ (* 2 i3)) 1.0)
			  sum (/ sign
				 (* k1 k2 k3 (+ (* k1 k1) (* k2 k2) (* k3 k3)))))))))))


(defun model-problem-computation (domain &key output plot)
  "~A - Solve the Laplace equation on a ~A.

Solves the Laplace problem with rhs identical 1 on the given domain.  The
solution strategy does uniform refinement and terminates if more than 20
seconds have passed after a step."
  (defparameter *result*
    (solve (blackboard
	    :problem (cdr-model-problem domain)
	    :output output :plot-mesh t
	    :success-if `(> :time 5.0))))
  (when plot
    (plot (getbb *result* :solution))))

(defun make-model-problem-demo (domain domain-name)
  (multiple-value-bind (title short long)
      (extract-demo-strings (documentation 'model-problem-computation 'function))
    (let ((demo
	   (make-demo
	    :name (format nil title domain-name)
	    :short (format nil short domain-name)
	    :long long
	    :execute (lambda ()
		       (model-problem-computation
			domain :plot t :output t)))))
      (adjoin-demo demo *laplace-demo*))))

(make-model-problem-demo (n-simplex-domain 1) "unit-interval")
(make-model-problem-demo (n-simplex-domain 2) "unit-triangle")
(make-model-problem-demo (n-cube-domain 2) "unit-quadrangle")
(make-model-problem-demo (n-simplex-domain 2) "unit-tetrahedron")
(make-model-problem-demo (tensorial-domain '(1 2)) "unit-wedge-1-2")
(make-model-problem-demo (tensorial-domain '(2 1)) "unit-wedge-2-1")
(make-model-problem-demo (n-cube-domain 3) "unit-cube")

(defun test-laplace-model-problem ()
  ;;; Elementary testing
  (let* ((dim 2) (order 1) (level 2)
	 (problem (cdr-model-problem dim))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (multiple-value-bind (matrix rhs)
	(discretize-globally problem h-mesh fedisc)
      (getrs (sparse-ldu matrix) rhs)))
  
  ;;; More testing
  (format t "~%~%*** Mesh convergence tests on laplace-test-problem ***~%")
  (let ((problem (cdr-model-problem 1)))
    (format t "~%1d-case (exact solution u(0.5) = 0.125)~%")
    (check-h-convergence problem 1 3 :order 1 :position #d(0.5))
    (check-h-convergence
     problem 1 3 :order 1 :position #d(0.5)
     :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t)
			    :success-if '(> :step 2))))
  (let ((problem (cdr-model-problem 2)))
    (format t "~%2d-case (exact solution u(0.5,0.5) = 0.0736713532...)~%")
    (check-h-convergence
     problem 1 4 :order 1 :position #d(0.5 0.5)
     :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t)
			    :success-if '(> :step 2)))
    (check-p-convergence problem 1 5 :level 1 :position #d(0.5 0.5)))
  (let ((problem (cdr-model-problem 3)))
    (format t "~%3d-case (exact solution u(0.5,0.5,0.5) = 0.0562128...)~%")
    (check-h-convergence
     problem 1 3 :order 1 :position #d(0.5 0.5 0.5)
     :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t)
			    :success-if '(> :step 2)))
    (check-p-convergence problem 1 4 :level 0 :position #d(0.5 0.5 0.5)))
  )

;;; (FL.application::test-laplace-model-problem)
(fl.tests:adjoin-test 'test-laplace-model-problem)
