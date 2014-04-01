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

;;; -*- mode: Lisp; mode : outline-minor; -*-

(in-package :fl.application)

;;; * Linear iterations

;;; ** Experiment 1.7

;;; *** Convergence of different methods on a small matrix
#+(or)
(let ((A #m((2.0 -1.0  0.0)
	    (-1.0  2.0 -1.0)
	    ( 0.0 -1.0  2.0)))
      (b #m((1.0) (2.0) (1.0))))
 ;(linsolve A b :output t :iteration *undamped-jacobi*)
  ;(linsolve A b :output t :iteration *gauss-seidel*)
  (linsolve A b :output t :iteration (make-instance '<sor> :omega 1.18))
  ;(linsolve A b :output t :iteration (make-instance '<gradient-method>))
  ;(linsolve A b :output t :iteration (make-instance '<cg>))
  ;(linsolve A b :output t :iteration (make-instance '<ilu>))
  )

;;; *** Testing linear iterations on the model problem

#+(or)
(progn
  (iteration-test (make-instance '<cg>) :size 20 :maxsteps 10)
  (iteration-test *undamped-jacobi* :size 5 :maxsteps 20 :output t)
  (iteration-test *gauss-seidel* :size 20 :maxsteps 200 :output t)
  (iteration-test *gauss-seidel* :dim 2 :size 5 :maxsteps 10 :output t)
  (iteration-test (make-instance '<ilu>) :dim 2 :size 5 :maxsteps 10 :output t)
  (iteration-test (make-instance '<cg>) :dim 2 :size 5 :maxsteps 10 :output t)
  (iteration-test (make-instance '<cg> :preconditioner *standard-ilu*)
		  :dim 2 :size 10 :output t :maxsteps 10))

(defun size-effort-graph (linit &key (dim 1) (from 3) (to 10) (maxsteps 100))
  "Returns a graph for the effort of linear iterations."
  (loop for size from from to to  ; :-)
	collecting
	(vector size (/ -1.0 (log (iteration-test linit :dim dim :size size :maxsteps maxsteps))))))

;;; *** Gauss-Seidel is two times faster than Jacobi
#+(or)
(gnuplot-polygons
 (list
  (cons "gs" (size-effort-graph *gauss-seidel* :maxsteps 50))
  (cons "jac" (size-effort-graph *undamped-jacobi* :maxsteps 50))
  ;;(cons "gradient" (size-effort-graph (make-instance '<gradient-method>)))
  ))

;;; *** Finding the optimal parameter for SOR

(defun sor-omega-graph (size &key (from 1.0) (to 1.9) (by 0.1))
  (loop for omega from from to to by by collect
	(vector omega (iteration-test (make-instance
				       '<sor> :omega omega) :size size))))
#+(or)
(gnuplot-polygon (sor-omega-graph 10 :from 1.5 :to 2.0 :by 0.1))

;;; Figure sor-cr
#+(or)
(gnuplot-polygons
 (list
  (cons "h=1/10" (sor-omega-graph 10 :from 1.5 :to 2.0 :by 0.01))
  (cons "h=1/20" (sor-omega-graph 20 :from 1.5 :to 2.0 :by 0.01)))
 :plot :file)

;;; * Multigrid

;;; ** Experiment 2.1

;;; *** Smoothing behaviour of simple smoothers

#+(or)
(multiple-value-bind (A b)
    #+(or)
    (model-problem-discretization :dim 2 :level 4)
    (model-problem-discretization :level 5)
  (let ((x (copy b)))
    (x<-0 b) (fill-random! x 1.0)
    (loop repeat 20 do
	  (plot x) (sleep 0.5)
	  (linsolve A b :sol x :output t :iteration
		    (ecase 'djac
		      (gs *gauss-seidel*)
		      (ujac *undamped-jacobi*)
		      (djac (make-instance '<jacobi> :damp 0.5))
		      )
		    :maxsteps 1))))

;;; ** Experiment: behaviour of the two-grid-algorithm

#+(or) ; try combinations dim/level/duration = 1/5/1.0 or 2/4/0.0
(let* ((dim 1) (level 5) (duration 1.0)
       (problem (cdr-model-problem dim))
       (cgc (geometric-cs :gamma 1 :pre-steps 0 :post-steps 0
			  :base-level (1- level))))
  (multiple-value-bind (A b)
      (problem-discretization problem :level level :order 1)
    (let ((x (copy b)))
      (x<-0 b) (fill-random! x 1.0)
      (loop repeat 10 do
	    (plot x) (sleep duration)
	    (linsolve A b :sol x :output t :iteration *gauss-seidel* :maxsteps 3)
	    (plot x) (sleep duration)
	    (linsolve A b :sol x :output t :iteration cgc :maxsteps 1)))))

;;; ** Experiment 2.19 (should be FMG, is multigrid is O(N))

#+(or)  ; dim=1: test with levels 1,2,3,..., dim=2: test with levels 1,2,3,4
(time
 (let* ((dim 2) (level 5)
	(problem (cdr-model-problem-on-domain dim))
	(geomg (geometric-cs :gamma 1 :pre-steps 0 :base-level 1)))
   (multiple-value-bind (A b)
       (problem-discretization problem :level level :order 1)
     (linsolve A b :iteration geomg :maxsteps 10))))

;;; ** Experiment: Inlay problem

#+(or)
(progn
  #+(or)
  (plot-diffusion (inlay-cell-problem 2 0.1) :refinements 2 :depth 0)
  ;; first order approximation
  (effective-tensor
   (cell-solve (inlay-cell-problem 2 0.1 0) :level 1
	       :order 1 :parametric (lagrange-mapping 2))))


;;; ** Multigrid convergence for the p-method with block smoothers

#+(or)  ; dim/level = 1/4, or 2/2: orders = 1,2,...
(let* ((dim 2) (level 2)
       (problem (cdr-model-problem dim))
       (geomg (geometric-cs :gamma 1 :pre-steps 0 :post-steps 2
			     :base-level 1)))
  (multiple-value-bind (A b)
      (problem-discretization problem :level level :order 1)
    (linsolve A b :output t :iteration geomg :maxsteps 20)))


;;; * Finite element methods

;;; ** convection problems
(defun convection-problem (dim eps)
  (cdr-model-problem dim
   :diffusion (scalar-diffusion dim eps)
   :convection (make-real-matrix (unit-vector dim 0))))

;;; *** dominant convection, choose eps = 1, 1/10, 1/64, 1/128, 1/200, 1/1000
#+(or)
(let ((dim 1))
  (plot (solve-laplace (convection-problem dim 1/1000) 4 1)))

;;; *** "stabilization" by using higher degree polynomials
#+(or)
(let* ((dim 1)
       (problem (convection-problem dim 1/1000))
       (sol (solve-laplace problem 5 2)))
  (plot sol)
  (loop for x from 0.9 upto 1.0 by 0.01 do
	(format t "x=~10,5F val = ~10,5F~%" x (fe-value sol (vector x))))
  )

;;; *** multigrid has problems only for large convection

#+(or)  ; dim=1: test with levels 1,2,3,..., dim=2: test with levels 1,2,3,4
(time
 (let* ((dim 1) (level 3)
	(problem (convection-problem dim 1/100))
	(geomg (geometric-cs :gamma 1 :pre-steps 0 :base-level 1)))
   (multiple-value-bind (A b)
       (problem-discretization problem :level level :order 1)
     (linsolve A b :iteration geomg :maxsteps 100))))




