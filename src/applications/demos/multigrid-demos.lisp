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

(in-package :fl.application)

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

(defun plot-iteration-behavior (dim order level iteration &key (nr-steps 10))
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

(defun make-plot-iteration-behavior-demo (dim order level iteration &key it-name)
  (unless it-name (setq it-name (class-name (class-of iteration))))
  (let ((demo
	 (make-demo
	  :name (format nil "~A-error-O~D-~DD" it-name order dim)
	  :short (format nil "Error development for ~A (order=~D, d=~D)." it-name order dim)
	  :long (format nil "~A~%Parameters: dim=~D, order=~D, level=~D, iteration=~A~%"
			(documentation 'plot-iteration-behavior 'function)
			dim order level it-name)
	  :execute
	  (lambda ()
	    (plot-iteration-behavior
	     dim order level
	     #'(lambda (A b x)
		 (linsolve A b :sol x :output t :iteration iteration :maxsteps 1)))))))
    (adjoin-demo demo *multigrid-demo*)))

(make-plot-iteration-behavior-demo 1 1 5 *gauss-seidel* :it-name "GS")
(make-plot-iteration-behavior-demo 2 1 3 *gauss-seidel* :it-name "GS")
(make-plot-iteration-behavior-demo 1 3 3 *gauss-seidel* :it-name "GS")

(defun make-two-grid-behavior-demo (dim order level)
  (let* ((cgc (geometric-cs :gamma 1 :pre-steps 0 :post-steps 0
			    :base-level (1- level)))
	 (smooth (make-instance '<multi-iteration> :base *gauss-seidel* :nr-steps 3))
	 (demo
	  (make-demo
	   :name (format nil "two-grid-method-~DD" dim)
	   :short (format nil "Error development for a two-grid method.")
	   :long (format nil "~A~%Parameters: dim=~D, order=~D, level=~D, smoother=GS~%"
			 (documentation 'plot-iteration-behavior 'function)
			 dim order level)
	   :execute
	   (lambda ()
	     (plot-iteration-behavior
	      dim order level
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

(defun smoother-performance-test (&key dim (order 1) (level 2) smoother
				  (output 1) simplex)
  "Tests performance of smoother on a Laplace model problem.
See make-smoother-demo for more information."
  (let* ((problem (cdr-model-problem (if simplex (n-simplex-domain dim) dim)))
	 (mm (uniformly-refined-hierarchical-mesh
	      (domain problem) level))
	 (fe-class (lagrange-fe order)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mm fe-class)
       (let ((result
	      (solve (blackboard :problem (lse :matrix mat :rhs rhs) :solver
				 (make-instance '<linear-solver> :iteration smoother
						:success-if '(> :step 10))
				 :output output))))
       (values (getbb result :step-reduction) result)))))

#+(or)
(let ((dim 2) (order 3))
  (smoother-performance-test :dim dim :order order :smoother (geometric-ssc)))

(defun smoother-demo-execute (smoother)
  (lambda ()
    (loop for dim = (user-input "Dimension (1..4): " #'parse-integer
				#'(lambda (x) (< 0 x 5)))
	  until (eq dim :up) do
	  (loop with max-order = (case dim (1 8) (2 7) (t 4))
		for order =
		(user-input (format nil "Order (1..~A): " max-order) #'parse-integer
			    #'(lambda (x) (<= 1 x max-order)))
		until (eq order :up) do
		(smoother-performance-test
		 :dim dim :order order :smoother smoother :output t)))))

(defun make-smoother-demo (smoother smoother-name)
  (let ((title (format nil "~A-performance" smoother-name))
	(short (format nil "Tests smoother ~A" smoother-name))
	(long (format nil "Tests the performance of ~A applied
to discretizations of different order of a Laplace model problem
on cubes of different dimensions." smoother-name)))
    (let ((demo
	   (make-demo
	    :name title :short short :long long
	    :execute (smoother-demo-execute smoother)
	    :test-input (format nil "2~%4~%up~%up~%"))))
      (adjoin-demo demo *multigrid-demo*))))

#+(or)(make-smoother-demo *gauss-seidel* "GS")
#+(or)(make-smoother-demo (geometric-ssc) "VC-SSC")

(defun smoother-graph-execute (smoother)
  (lambda ()
    (loop for dim = (user-input "Dimension (1..4): " #'parse-integer
				#'(lambda (x) (< 0 x 5)))
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
  (let ((title (format nil "~A-cr-graph" smoother-name))
	(short (format nil "Plots graph 'order->CR(order)' for ~A" smoother-name))
	(long (format nil "Plots a graph of the convergence rate
for ~A smoother applied to discretizations of different order of
a Laplace model problem on cubes of different dimensions." smoother-name)))
    (let ((demo
	   (make-demo
	    :name title :short short :long long
	    :execute (smoother-graph-execute smoother)
	    :test-input (format nil "2~%up~%"))))
      (adjoin-demo demo *multigrid-demo*))))

(make-smoother-performance-graph-demo *gauss-seidel* "GS")
(make-smoother-performance-graph-demo (geometric-ssc) "VC-SSC")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; BPX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun bpx-demo-computation (problem order level &optional (smoother-class '<jacobi>))
  "Performs the BPX demo."
  (let* ((smoother (make-instance smoother-class :damp 1.0))
	 (cs (geometric-cs
	      :gamma 1 :smoother smoother :pre-steps 1 :post-steps 0
	      :combination :additive :base-level 1))
	 (bpx (make-instance '<cg> :preconditioner cs))
	 (solver (make-instance '<linear-solver> :iteration bpx
				:success-if `(or (< :defnorm 1.0e-12) (> :step 20))))
	 (mesh (uniformly-refined-hierarchical-mesh (domain problem) level))
	 (nr-comps (nr-of-components problem))
	 (fedisc (lagrange-fe order :nr-comps nr-comps)))
    (multiple-value-bind (A b)
	(discretize-globally problem mesh fedisc)
      (setq *result*
	    (solve (blackboard :matrix A :rhs b :solver solver :output 1))))))

;;; (bpx-demo-computation (cdr-model-problem 2) 1 2 '<jacobi>)
;;; (bpx-demo-computation (cdr-model-problem 2) 1 2 '<psc>)
;;; (bpx-demo-computation (elasticity-model-problem 2) 5 2 '<jacobi>)
;;; (bpx-demo-computation (elasticity-model-problem 2) 5 2 '<geometric-psc>)

(defun make-bpx-demo (problem problem-name order level)
  (adjoin-demo
   (make-demo
    :name (format nil "BPX-~A" problem-name)
    :short (format nil "BPX for the problem ~A" problem-name)
    :long (format nil "Shows convergence for BPX for the
problem ~A.~%Parameters: dim=~D, order=~D, level=~D,~%"
		  problem-name (dimension (domain problem)) order level)
    :execute (lambda () (bpx-demo-computation problem order level)))
   *multigrid-demo*))
(make-bpx-demo (cdr-model-problem 2) "laplace-on-square" 1 4)
(make-bpx-demo (elasticity-model-problem 2) "elasticity-on-square" 1 4)


;;;; Testing:

(defun test-multigrid-demos ()
  (smoother-performance-test :dim 1 :order 1 :smoother *gauss-seidel* :output t)
  (smoother-performance-test :dim 1 :order 1 :smoother (geometric-psc))
  (smoother-performance-test :dim 1 :order 6 :smoother (geometric-ssc))
  (smoother-performance-test :dim 3 :order 4 :level 2 :simplex t :smoother (geometric-ssc))
  (bpx-demo-computation (cdr-model-problem 2) 1 5 '<jacobi>)
  (time (bpx-demo-computation (elasticity-model-problem 3) 5 2))
  )

;;; (test-multigrid-demos)
(fl.tests:adjoin-test 'test-multigrid-demos)
