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

(defun further-laplace-tests ()
  "This function provides further tests for Laplace problems, partially on
non-structured meshes and/or on domains with curved boundary."

  ;; 1D tests
  (format t "~%1d Laplace test with nonzero rhs~%")
  
  (let ((problem
	 (let ((dim 1))
	   (cdr-model-problem (n-cube-domain dim)
			      :source #'(lambda (x) #I"exp(x[0])")))))
    (check-h-convergence problem 1 6 :order 1 :position #d(0.5))
    (check-p-convergence problem 1 10 :level 1 :position #d(0.5)))
  
  (format t "~%1d Laplace test with Dirichlet bc~%")
  (let* ((dim 1)
	 (problem
	  (cdr-model-problem
	   (n-cube-domain dim)
	   :source #'(lambda (x) #I"exp(x[0])")
	   :dirichlet #'(lambda (x)
			  (if (zerop (aref x 0)) 0.0 1.0)))))
    (check-h-convergence problem dim 6 :order 1 :position #d(0.5))
    (check-p-convergence problem dim 5 :level 0 :position #d(0.5)))

  ;; 2D tests
  (format t "~%2d Laplace test on a triangle domain~%")
  
  (let ((problem (cdr-model-problem (n-simplex-domain 2))))
    (check-h-convergence
     problem 2 6 :order 1 :position #d(0.33 0.33)
     :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t)
			    :success-if '(> :step 2)))
    (check-p-convergence problem 2 6 :level 2 :position #d(0.33 0.33)))
  
  (format t "~%2d Laplace test on the unit circle, exact solution u(0,0) = 0.25~%")
  (let ((problem (cdr-model-problem (n-ball-domain 2))))
    (format t "~%h-refinement, O(h^2) convergence~%")
    (check-h-convergence
     problem 1 3 :order 1 :position #d(0.0 0.0)
     :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t)
			    :success-if '(> :step 2)))
    (format t "~%p-refinement, no convergence, because domain is not approximated~%")
    (check-p-convergence problem 1 3 :level 2 :position #d(0.0 0.0))
    (format t "~%p-refinement, exponential convergence due to isoparametric approximation~%~
Exact value: u(0.1,0.1)=0.245~%")
    ;; We do not use the origin as ealuation point in the following line
    ;; because it sometimes is not recognized as being contained in any
    ;; cell due to approximation errors when solving the nonlinear mappings
    ;; (interestingly, I observed this only with Allegro CL).  Maybe a
    ;; search in lower-dimensional cells could help?
    (check-p-convergence problem 1 5 :level 2 :position #d(0.1 0.1) :isopar t))

  (format t "~%Laplace with exact solution u=exp(x+y), -> u(1/2,1/2)=1.6487212707~%")
  (let* ((dim 2)
	 (problem
	  (cdr-model-problem
	   (n-cube-domain dim)
	   :source #'(lambda (x) #I"-2.0*exp(x[0]+x[1])")
	   :dirichlet #'(lambda (x) #I"exp(x[0]+x[1])"))))
    (check-h-convergence
     problem 1 3 :order 1 :position #d(0.25 0.25)
     :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t)
			    :success-if '(> :step 2)))
    (check-p-convergence problem 1 5 :level 0 :position #d(0.25 0.25)))
  
  ;; 3D tests
  (format t "~%3d Laplace test on a tetrahedron~%")
  (time
   (let ((problem (cdr-model-problem (n-simplex-domain 3))))
     (check-h-convergence
      problem 2 3 :order 1 :position #d(0.25 0.25 0.25)
     :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t :base-level 2)
			    :success-if '(> :step 2)))
     (check-p-convergence problem 1 3 :level 2 :position #d(0.25 0.25 0.25))))

  (format t "~%Laplace with exact solution u=x*y*z(1-x-y-z), i.e. u(1/4,1/4,1/4)=1/256=3.90625e-3~%")
  (let* ((domain (n-simplex-domain 3))
	 (problem
	  (cdr-model-problem
	   domain :source #'(lambda (x) #I(2.0*(x[1]*x[2]+x[0]*x[2]+x[0]*x[1]))))))
    (check-h-convergence
     problem 2 4 :order 1 :position #d(0.25 0.25 0.25)
     :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t :base-level 2)
			    :success-if '(> :step 2)))
    ;; we should obtain the exact solution for ansatz spaces that include
    ;; the solution
    (check-p-convergence problem 1 5 :level 0 :position #d(0.25 0.25 0.25))
    (check-p-convergence problem 4 4 :level 1 :position #d(0.25 0.25 0.25)))

  (format t "~%3d Laplace test on the unit ball, exact solution u(0,0,0)=1/6=0.1666...~%")
  (let ((problem (cdr-model-problem (n-ball-domain 3))))
    (time (check-h-convergence
	   problem 0 3 :order 1 :position #d(0.0 0.0 0.0)
	   :solver (make-instance '<linear-solver> :iteration (geometric-cs :fmg t)
				  :success-if '(> :step 2))))
    )
  )

;;; (fl.application::further-laplace-tests)
(fl.tests:adjoin-test 'further-laplace-tests)

