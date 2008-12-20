;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; evp-cdr.lisp - computing eigenvalues for cdr problems
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

(file-documentation "This file contains routines for computing
eigenvalue/eigenvector pairs for convection-diffusion-reaction problems.")

(defvar *laplace-eigenvalue-demo*
  (make-demo
   :name "eigenvalues"
   :short "Some eigenvalue problems"))
(adjoin-demo *laplace-eigenvalue-demo* *laplace-demo*)

(defun laplace-eigenvalue-computation
    (domain &key output plot (dirichlet 0.0))
  "Function performing the eigenvalue demo for the Laplace operator."
  (let ((problem (cdr-model-problem
		  domain :evp t :dirichlet dirichlet)))
    (storing
      (solve (blackboard
	      :problem problem :base-level 3
	      :success-if '(or (>= :time 20) (>= :nr-levels 5))
	      :output output :observe
	      (append *stationary-fe-strategy-observe*
		      (list (list "             lambda" "~19,10,2E"
				  #'(lambda (bb) (declare (ignore bb))
					    (unbox (slot-value problem 'eigenvalues)))))))))
    (when plot
      (plot (getbb *result* :solution)))))

(defun make-laplace-eigenvalue-demo (domain domain-name)
  (let ((title domain-name)
	(short (format nil "Eigenvalues of Laplace on a ~A." domain-name))
	(long (format nil "Computes eigenvalues for the Laplace
operator on a ~A.  The solution strategy does uniform refinement
and terminates if more than 20 seconds have passed after a
step." domain-name)))
    (let ((demo
	   (make-demo
	    :name title :short short :long long
	    :execute (lambda ()
		       (laplace-eigenvalue-computation
			domain :plot t :output 1)))))
      (adjoin-demo demo *laplace-eigenvalue-demo*))))

(make-laplace-eigenvalue-demo (n-simplex-domain 1) "unit-interval")
(make-laplace-eigenvalue-demo (n-cube-domain 2) "unit-quadrangle")

#+(or)
(let ((problem (cdr-model-problem
		(n-cube-domain 1)
		:multiplicity 2
		:reaction #m(1.0) :source nil :dirichlet nil :evp t)))
  (storing
    (solve (blackboard
	    :problem problem :base-level 2
	    :success-if '(>= :nr-levels 4)
	    :output :all :observe
	    (append *stationary-fe-strategy-observe*
		    (list (list "             lambda" "~A"
				#'(lambda (bb)
				    (declare (ignore bb))
				    (slot-value problem 'eigenvalues)))))))))

;;;; Testing

(defun evp-cdr-test ()
  (plot (laplace-eigenvalue-computation
	 (n-cube-domain 1) :output t :plot t))

  ;; the following is used in the manual
  (let ((problem (cdr-model-problem 2 :evp t)))
    (storing
      (solve (blackboard :problem problem
			 :success-if '(or (>= :time 5) (>= :nr-levels 5))
			 :output 1))))
  (slot-value (getbb *result* :problem) 'eigenvalues)
  (plot (getbb *result* :solution))
  
  ;; bug shows up when using Lispworks
  (let ((A #m((8.0 -12.0) (-12.0 18.0)))
        (B (eye 2)))
    (hegv A B :V))
  
  )

;;; (evp-cdr-test)
(fl.tests:adjoin-test 'evp-cdr-test)