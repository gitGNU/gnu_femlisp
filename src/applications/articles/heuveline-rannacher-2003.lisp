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

(in-package :fl.application)

(defparameter *HR-evaluation-point*
  #d(0.5 2.5))

(defparameter *HR-gradx-value*
  (* pi 0.5 (sin (* 0.75 pi)) (sin (* 2.625 pi))))

(defun heuveline-rannacher-rhs (x)
  "Smooth rhs for the problem -Delta u = rhs in the H-R article."
  #I(-0.8125*pi*pi * sin(0.5*pi*(x[0]+1.0)) * sin(0.75*pi*(x[1]+1.0))))

(defun heuveline-rannacher-domain ()
  (box-domain '((-1.0 1.0) (-1.0 3.0))))

(defun heuveline-rannacher-problem ()
  (let ((domain (heuveline-rannacher-domain)))
    (cdr-model-problem
     domain :source #'heuveline-rannacher-rhs)))

#+(or)
(check (uniform-mesh-on-box-domain (heuveline-rannacher-domain) #(1 2)))

(defparameter *HR-delta* 1.0e-12
  "Small positive value which is used for ``smoothing'' the distributional
rhs.")

(defun heuveline-rannacher-dual-problem-rhs (&key cell fe &allow-other-keys)
  "Distributional rhs for the dual problem of the article by
Heuveline&Rannacher.  We distribute it to several points to
ensure that all surrounding cells contribute.  Warning: This
function assumes that a structured cube mesh is used!"
  (let ((rhs (make-local-vec fe))
	(local (global->local cell *HR-evaluation-point*))
	(delta *HR-delta*))
    (when (every #'(lambda (c) (<= (- delta) c (+ 1.0 delta))) local)
      ;; point is inside or very near the cell...
      (let ((nr-neighbors 
	     (expt 2 (count-if #'(lambda (c)
				   (or (< (abs c) delta)
				       (< (abs (1- c)) delta)))
			       local)))
	    (Dphi^-1 (m/ (local->Dglobal cell local))))
	(loop+ (i (shape (fe-basis fe))) do
	   (setf (vref rhs i)
		 (/ (vref (m* (ensure-matlisp (coerce (evaluate-gradient shape local) 'double-vec)
					      :row)
			      Dphi^-1)
			  0)
		    nr-neighbors)))))
    rhs))

(defun heuveline-rannacher-dual-problem-fe-rhs ()
  #'(lambda (cell)
      (when (= (dimension cell) 2)
	(list 'FL.CDR::FE-RHS
	      (make-instance
	       '<coefficient> :demands '(:cell :fe) :eval
	       #'heuveline-rannacher-dual-problem-rhs)))))

(defun heuveline-rannacher-dual-problem ()
  "Only for testing purposes, otherwise such a problem is
automatically generated inside the error estimator cycle."
  (dual-problem (heuveline-rannacher-problem)
		(heuveline-rannacher-dual-problem-fe-rhs)))

(defun heuveline-rannacher-computation (order levels &key (output 1) plot)
  "Performs the Heuveline-Rannacher demo."
  (defparameter *result*
    (let ((*output-depth* output))
      (solve
       (make-instance
	'<stationary-fe-strategy>
	:fe-class (lagrange-fe order)
	:estimator (make-instance '<duality-error-estimator> :functional
				  (heuveline-rannacher-dual-problem-fe-rhs))
	:indicator (make-instance '<largest-eta-indicator>
				  :pivot-factor 0.2 :from-level 1)
	:success-if `(or (<= :global-eta 1.0e-10) (= :max-level ,(1- levels)))
	:plot-mesh plot :observe
	(append *stationary-fe-strategy-observe*
		(list
		 (list " grad-x (0.5,2.5)" "~17,10,2E"
		       #'(lambda (blackboard)
			   (with-items (&key solution) blackboard
			     (vref (aref (fe-gradient solution *HR-evaluation-point*) 0) 0))))
		 *eta-observe*)))
       (let ((problem (heuveline-rannacher-problem)))
	 (blackboard :problem problem :mesh
		     (change-class (uniform-mesh-on-box-domain (domain problem) #(1 2))
				   '<hierarchical-mesh>))))))
    (when plot
      (plot (getbb *result* :solution) :depth 3)))

#+(or) (heuveline-rannacher-computation 4 4 :output :all :plot t)

(defun make-heuveline-rannacher-demo (order levels)
  (let ((title "HR2003-1")
	(short "Solves problem 1 in [Heuveline-Rannacher 2003]")
	(long (format nil "Solves for the functional
@math{J(phi)=dphi/dx(0.5,2.5)} evaluated for approximations
@math{u_k} to a Poisson equation with exact solution:

@math{u(x,y) = sin(pi/2*(x+1))*sin(3/4*pi*(y+1))}

Thus, the precise value is 1.026172152977031.
Parameters of the computation: order=~D, levels=~D." order levels)))
    (let ((demo
	   (make-demo
	    :name title :short short :long long
	    :execute (lambda ()
		       (heuveline-rannacher-computation
			order levels :output 1 :plot t)))))
      (adjoin-demo demo *articles-demo*))))

(make-heuveline-rannacher-demo 4 5)

#|
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
       (defparameter *result* (getrs (sparse-ldu mat) rhs)))))
  (plot *result*)

  (let ((gradx (vref (vref (fe-gradient *result* *HR-evaluation-point*) 0) 0)))
    (list gradx *HR-gradx-value* (- gradx *HR-gradx-value*)))

  (time
   (let* ((order 4) (level 3)
	  (problem (heuveline-rannacher-dual-problem))
	  (fe-class (lagrange-fe order))
	  (mesh (uniformly-refined-hierarchical-mesh (domain problem) level)))
     (multiple-value-bind (mat rhs)
	 (discretize-globally problem mesh fe-class)
       (defparameter *result* (getrs (sparse-ldu mat) rhs)))))
  (plot *result*)

  )  
