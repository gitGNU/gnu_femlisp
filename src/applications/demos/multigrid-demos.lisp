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
  (let* ((problem (cdr-model-problem dim))
	 (h-mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (fedisc (lagrange-fe order)))
    (multiple-value-bind (A b)
	(discretize-globally problem h-mesh fedisc)
      (let ((x (copy b)))
	(x<-0 b) (fill-random! x 1.0)
	(loop repeat nr-steps do
	      (plot x) (sleep 0.5)
	      (funcall iteration A b x)
	      finally (plot x))))))

(defun make-plot-iteration-behavior-demo (dim level order iteration &key it-name)
  (unless it-name (setq it-name (class-name (class-of iteration))))
  (let ((demo
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

(make-plot-iteration-behavior-demo 1 5 1 *gauss-seidel* :it-name "GS")

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; Robust smoothing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun smoother-performance-test (&key dim order (level 2) smoother
				  output simplex)
  "Tests performance of smoother on a Laplace model problem.
See make-smoother-demo for more information."
  (let* ((problem (cdr-model-problem (if simplex (n-simplex-domain dim) dim)))
	 (mm (uniformly-refined-hierarchical-mesh
	      (domain problem) level))
	 (fe-class (lagrange-fe order)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mm fe-class)
       (let ((result
	      (solve (make-instance '<linear-solver> :iteration smoother
				    :success-if '(> :step 10) :output output)
		     (blackboard :matrix mat :rhs rhs))))
       (values (getbb result :step-reduction) result)))))

#+(or)
(let ((dim 1) (order 6))
  (smoother-performance-test :dim dim :order order :smoother *gauss-seidel*))

(defun smoother-demo-execute (smoother)
  (lambda ()
    (loop for dim = (user-input "Dimension (1..4): "
				#'(lambda (x) (and (integerp x) (< 0 x 5))))
	  until (eq dim :up) do
	  (loop with max-order = (case dim (1 8) (2 7) (t 4))
		for order =
		(user-input (format nil "Order (1..~A): " max-order)
			    #'(lambda (x) (and (integerp x) (<= 1 x max-order))))
		until (eq order :up) do
		(smoother-performance-test
		 :dim dim :order order :smoother smoother :output t)))))

(defun make-smoother-demo (smoother smoother-name)
  "~name~-performance - Tests smoother ~name~

Tests the performance of ~name~ applied to discretizations of
different order of a Laplace model problem on cubes of different
dimensions."
  (multiple-value-bind (name short long)
      (extract-demo-strings
       (documentation 'make-smoother-demo 'function)
       (list (cons "~name~" smoother-name)))
    (let ((demo
	   (make-demo
	    :name name :short short :long long
	    :execute (smoother-demo-execute smoother)
	    :test-input (format nil "2~%4~%up~%up~%"))))
      (adjoin-demo demo *multigrid-demo*))))

#+(or)(make-smoother-demo *gauss-seidel* "GS")
#+(or)(make-smoother-demo (geometric-ssc) "VC-SSC")

(defun smoother-graph-execute (smoother)
  (lambda ()
    (loop for dim = (user-input "Dimension (1..4): "
				#'(lambda (x) (and (integerp x) (< 0 x 5))))
	  until (eq dim :up) do
	  (plot
	   (list
	    (cons
	     "convergence-rate"
	     (loop for order from 1 upto (case dim (1 8) (2 4) (t 3))
		   collect
		   (vector order
			   (smoother-performance-test
			    :dim dim :order order :smoother smoother)))))
	   :left 0 :right 8 :top 1.0 :bottom 0.0
	   ))))

(defun make-smoother-performance-graph-demo (smoother smoother-name)
  "~name~-cr-graph - Plots graph 'order->CR(order)' for ~name~

Plots a graph of the convergence rate for ~name~ smoother
applied to discretizations of different order of a Laplace model
problem on cubes of different dimensions."
  (multiple-value-bind (name short long)
      (extract-demo-strings
       (documentation 'make-smoother-performance-graph-demo 'function)
       (list (cons "~name~" smoother-name)))
    (let ((demo
	   (make-demo
	    :name name :short short :long long
	    :execute (smoother-graph-execute smoother)
	    :test-input (format nil "2~%up~%"))))
      (adjoin-demo demo *multigrid-demo*))))

(make-smoother-performance-graph-demo *gauss-seidel* "GS")
(make-smoother-performance-graph-demo (geometric-ssc) "VC-SSC")

(defun test-multigrid-demos ()
  (smoother-performance-test :dim 1 :order 6 :smoother (geometric-ssc))
  (smoother-performance-test :dim 3 :order 4 :level 3 :simplex t :smoother (geometric-ssc))
  )
