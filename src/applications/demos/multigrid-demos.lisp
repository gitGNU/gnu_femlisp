;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; multigrid-demos.lisp
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

(defvar *multigrid-demo*
  (make-demo
   :name "multigrid"
   :short "Demos focused on multigrid as solution technique"
   :long
   "Multigrid is a fast solution technique for problems with
elliptic character.  Here are some demos showing the application
of multigrid."))

(adjoin-demo *multigrid-demo* *solver-demo*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Smoothing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun plot-iteration-behavior (dim level order iteration &key (nr-steps 10))
  "Plots the propagation of the error for some iteration.  The
iteration is applied to a finite element discretization with
Lagrange finite elements of a given order of the Laplace model
problem in dim space dimensions on a mesh with mesh-width
2^{-level}.  This is done by plotting the approximation starting
with an initial random guess and right-hand side 0."
  (let* ((problem (laplace-test-problem-on-domain (n-cube-domain dim)))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (multiple-value-bind (A b)
	(discretize-globally problem h-mesh fedisc)
      (let ((x (copy b)))
	(x<-0 b) (x<-random x 1.0d0)
	(loop repeat nr-steps do
	      (plot x) (sleep 0.5)
	      (funcall iteration A b x)
	      finally (plot x))))))

(defun make-plot-iteration-behavior-demo (dim level order iteration)
  (let* ((it-name (class-name (class-of iteration)))
	 (demo
	  (make-demo
	   :name (format nil "~A-error-~DD" it-name dim)
	   :short (format nil "Error development for ~A." it-name)
	   :long (format nil "~A~%Parameters: dim=~D, level=~D, order=~D, iteration=~A~%"
			 (documentation 'plot-iteration-behavior 'function)
			 dim level order it-name)
	   :execute
	   (lambda ()
	     (plot-iteration-behavior
	      dim level order
	      #'(lambda (A b x)
		  (linsolve A b :sol x :output t :iteration iteration :maxsteps 1)))))))
    (adjoin-demo demo *multigrid-demo*)))

(make-plot-iteration-behavior-demo 1 5 1 *gauss-seidel*)

(defun make-two-grid-behavior-demo (dim level order)
  (let* ((cgc (geometric-cs :gamma 1 :pre-steps 0 :post-steps 0
			    :base-level (1- level)))
	 (smooth (make-instance '<multi-iteration> :base *gauss-seidel* :nr-steps 3))
	 (demo
	  (make-demo
	   :name (format nil "two-grid-method-~DD" dim)
	   :short (format nil "Error development for a two-grid method.")
	   :long (format nil "~A~%Parameters: dim=~D, level=~D, order=~D, smoother=GS~%"
			 (documentation 'plot-iteration-behavior 'function)
			 dim level order)
	   :execute
	   (lambda ()
	     (plot-iteration-behavior
	      dim level order
	      (let ((cgc-p nil))
		#'(lambda (A b x)
		    (linsolve A b :sol x :output t :maxsteps 1
			      :iteration (if cgc-p cgc smooth))
		    (setq cgc-p (not cgc-p))))
	      :nr-steps 6)))))
    (adjoin-demo demo *multigrid-demo*)))

(make-two-grid-behavior-demo 1 5 1)
