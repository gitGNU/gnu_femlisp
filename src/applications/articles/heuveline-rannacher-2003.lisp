;;; -*- mode: lisp; fill-column: 64; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; heuveline-rannacher-2003.lisp
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

(defparameter *HR-evaluation-point*
  (double-vec 0.5d0 2.5d0))

(defparameter *HR-gradx-value*
  (* pi 0.5d0 (sin (* 0.75d0 pi)) (sin (* 2.625d0 pi))))

(defun heuveline-rannacher-rhs (x)
  "Smooth rhs for the problem -Delta u = rhs in the H-R article."
  #I(-0.8125d0*pi*pi * sin(0.5d0*pi*(x[0]+1.0d0)) * sin(0.75d0*pi*(x[1]+1.0d0))))

(defun heuveline-rannacher-domain ()
  (box-domain '((-1.0 1.0) (-1.0 3.0))))

(defun heuveline-rannacher-problem ()
  (let ((domain (heuveline-rannacher-domain)))
    (cdr-model-problem
     domain :source #'heuveline-rannacher-rhs)))

#+(or)
(check (uniform-mesh-on-box-domain (heuveline-rannacher-domain) #(1 2)))

(defparameter *HR-delta* 1.0d-12
  "Small positive value which is used for ``smoothing'' the distributional
rhs.")

(defun heuveline-rannacher-dual-problem-rhs (cell fe)
  "Distributional rhs for the dual problem of the
Heuveline-Rannacher article.  We distribute it to several points
to ensure that all surrounding cells contribute.  Warning: This
function assumes that a structured cube mesh is used!"
  (let ((rhs (make-local-vec fe))
	(local (global->local cell *HR-evaluation-point*))
	(delta *HR-delta*))
    (when (every #'(lambda (c) (<= (- delta) c (+ 1.0d0 delta))) local)
      ;; point is inside or very near the cell...
      (loop
       with nr-neighbors =
       (expt 2 (count-if #'(lambda (c)
			     (or (< (abs c) delta)
				 (< (abs (1- c)) delta)))
			 local))
       with Dphi^-1 = (m/ (local->Dglobal cell local))
       for shape in (fe-basis fe)
       and i from 0 do
       (setf (matrix-ref rhs i)
	     (/ (matrix-ref (m* (ensure-matlisp (evaluate-gradient shape local) :row)
				Dphi^-1) 0)
		nr-neighbors))))
    rhs))

(defun heuveline-rannacher-dual-problem-fe-rhs ()
  #'(lambda (cell)
      (when (= (dimension cell) 2)
	(list 'CDR::FE-RHS #'heuveline-rannacher-dual-problem-rhs))))
  
(defun heuveline-rannacher-dual-problem ()
  (dual-problem (heuveline-rannacher-problem)
		(heuveline-rannacher-dual-problem-fe-rhs)))

(defparameter *result* nil)

(defun heuveline-rannacher-computation (order levels &key output plot)
  "HR2003-1 - Solves problem 1 in [Heuveline-Rannacher 2003]

Solves for the functional $J(phi)=dphi/dx(0.5,2.5)$ evaluated
for approximations $u_k$ to a Poisson equation with exact
solution:

$$  u(x,y) = sin(pi/2*(x+1))*sin(3/4*pi*(y+1))  $$

Thus, the precise value is 1.026172152977031d0.
Parameters of the computation: order=~order~, levels=~levels~."
  (defparameter *result*
    (let ((problem (heuveline-rannacher-problem))
	  (solver
	   (?2 *lu-solver*
	       (make-instance
		'<linear-solver>
		:iteration
		(let ((smoother (geometric-ssc)))
		  (make-instance '<s1-reduction> :max-depth 2 :pre-steps 1 :pre-smooth smoother
				 :post-steps 1 :post-smooth smoother
				 :gamma 1 :coarse-grid-iteration
				 (?2 *lu-iteration*
				     (make-instance '<s1-coarse-grid-iterator>))))
		:success-if `(and (< :defnorm 1.0e-10) (> :step-reduction 0.9))
		:failure-if `(and (> :step 2) (> :step-reduction 0.9))
		:output (eq output :all)))))
      (solve
       (make-instance
	'<stationary-fe-strategy>
	:fe-class (lagrange-fe order) :solver solver
	:estimator (make-instance '<duality-error-estimator> :functional
				  (heuveline-rannacher-dual-problem-fe-rhs))
	:indicator (make-instance '<largest-eta-indicator>
				  :pivot-factor 0.2 :from-level 1)
	:success-if `(or (<= :global-eta 1.0e-10) (= :max-level ,(1- levels)))
	:plot-mesh plot	:output output :observe
	(append *stationary-fe-strategy-observe*
		(list
		 (list "   grad-x (0.5,2.5)" "~19,10,2E"
		       #'(lambda (blackboard)
			   (with-items (&key solution) blackboard
			     (vec-ref (aref (fe-gradient solution *HR-evaluation-point*) 0) 0))))
		 (list "         ETA" "~12,2,2E"
		       #'(lambda (blackboard)
			   (getbb blackboard :global-eta))))))
       (blackboard :problem problem
		   :mesh (change-class (uniform-mesh-on-box-domain (domain problem) #(1 2))
				       '<hierarchical-mesh>)
		   :output t))))
    (when plot
      (plot (getf *result* :solution) :depth 3)))



#+(or) (heuveline-rannacher-computation 4 6 :output t :plot t)

(defun make-heuveline-rannacher-demo (order levels)
  (multiple-value-bind (title short long)
      (extract-demo-strings
       (documentation 'heuveline-rannacher-computation 'function)
       `(("~order~" . ,order) ("~levels~" . ,levels)))
    (let ((demo
	   (make-demo
	    :name title :short short :long long
	    :execute (lambda ()
		       (heuveline-rannacher-computation
			order levels :output t :plot t)))))
      (adjoin-demo demo *articles-demo*))))

(make-heuveline-rannacher-demo 4 4)

#| Konvergenzraten durchweg gut...
* (heuveline-rannacher-computation 4 8 :output t)
auf toba?:
 CELLS     DOFS  MENTRIES   TIME   grad-x (0.5,2.5)                ETA
    2        45      1225    2.1   1.0098237387d+00   2.2392139252d-02
    8       198      4753   10.8   1.0246749201d+00   1.3097813462d-03
   11       279      7266   24.0   1.0260207607d+00   1.8238160389d-04
   23       568     15010   58.4   1.0261454565d+00   2.0900486250d-05
   32       793     23531  106.2   1.0261860710d+00   2.2922568999d-06
   50      1201     35959  192.6   1.0261714327d+00   8.5627095453d-07
   77      1834     57159  322.9   1.0261721980d+00   3.8457544838d-08
danach schlechtere Konvergenz und schlechtere Werte

auf ortler:
 CELLS     DOFS  MENTRIES   TIME   grad-x (0.5,2.5)                ETA
    2        45      1225    0.6   1.0098237387d+00   2.2392139252d-02
    8       198      4753    6.1   1.0246749201d+00   1.3097813462d-03
   11       279      7266   10.4   1.0260207607d+00   1.8238160389d-04
   23       568     15010   21.7   1.0261454566d+00   2.0900484053d-05
   32       793     23531   37.9   1.0261860710d+00   2.2922561027d-06
   50      1201     35959   67.8   1.0261714327d+00   8.5627035398d-07
   77      1834     57159  113.0   1.0261721980d+00   3.8444739754d-08
|#

(defun test-heuveline-rannacher ()

  (time
   (let* ((order 4) (level 2)
	  (problem (heuveline-rannacher-problem))
	  (fe-class (lagrange-fe order))
	  (mesh (uniformly-refined-hierarchical-mesh (domain problem) level)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mesh fe-class)
       (defparameter *result* (m* (sparse-ldu mat) rhs)))))
  (plot *result*)

  (let ((gradx (matrix-ref (fe-gradient *result* *HR-evaluation-point*) 0)))
    (list gradx *HR-gradx-value* (- gradx *HR-gradx-value*)))
  (fe-gradient *result* *HR-evaluation-point*)

  (time
   (let* ((order 4) (level 3)
	  (problem (heuveline-rannacher-dual-problem))
	  (fe-class (lagrange-fe order))
	  (mesh (uniformly-refined-hierarchical-mesh (domain problem) level)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mesh fe-class)
       (defparameter *result* (m* (sparse-ldu mat) rhs)))))
  (plot *result*)

  )  
  
  

