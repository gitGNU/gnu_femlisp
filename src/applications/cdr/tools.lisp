;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; tools.lisp - some tools (Warning: old code)
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

(defun problem-discretization (problem &key parametric (level 1) (order 1) &allow-other-keys)
  "Returns stiffness matrix and rhs for a discretization of order
@var{order} on a mesh with a mesh-size 2^@var{-level}."
  (let ((mm (uniformly-refined-hierarchical-mesh (domain problem) level
						  :parametric parametric))
	(fedisc (lagrange-fe order)))
    (discretize-globally problem mm fedisc)))

(defun model-problem-discretization (&key (dim 1) (size 1) (level 1) (order 1) &allow-other-keys)
  "Returns stiffness matrix and right-hand side for the discretization of
the Laplace model problem on a cube of dimension @var{dim} using a mesh of
meshwidth 1/@var{size} cells in each dimension with Lagrange finite
elements of order @var{order}."
  (let* ((problem (cdr-model-problem dim))
	 (fedisc (lagrange-fe order))
	 (domain (n-cube-domain dim))
	 (mm (change-class (uniform-mesh-on-box-domain domain size) '<hierarchical-mesh>)))
    (loop repeat level do (refine mm))
    (discretize-globally problem mm fedisc)))

(defun iteration-test (linit &rest args &key (maxsteps 200) output &allow-other-keys)
  "Tests the linear iteration @var{linit} on a model problem specified via
keyword parameters in @var{args}."
  (multiple-value-bind (A b constraints-P constraints-Q constraints-r)
      (apply #'model-problem-discretization args)
    (declare (ignore constraints-Q))
    (let ((x (copy b)))
      (x<-0 b) (fill-random! x 1.0)
      (fill-random! x 1.0)
      (x<-Ay x constraints-P constraints-r)
      (let ((ls (make-instance '<linear-solver> :iteration linit
			       :success-if `(or (< :defnorm 1.0e-12) (> :step ,maxsteps))
			       :output output)))
	(setq *result* (solve ls (blackboard :problem (lse :matrix A :rhs b :solution x))))
	(getbb *result* :report)))))

(defun solve-laplace (problem level order &key parametric (solver *lu-solver*))
  "An old routine for solving the Laplace problem.  You should use the
interface provided by @code{solve} which provides automatic and more
general solving together with a more flexible customization."
  (multiple-value-bind (matrix rhs)
      (problem-discretization problem :level level :order order :parametric parametric)
    (getbb (solve solver (blackboard :problem (lse :matrix matrix :rhs rhs)))
	   :solution)))

(defun check-h-convergence (problem min-level max-level
			    &key order position (solver *lu-solver*))
  (format t "~%Dimension = ~D, order = ~D, midpoint = ~A:~%"
	  (dimension (domain problem)) order position)
  (format t "Level :  u_h(midpoint)  :  (u_h-u_2h)(midpoint)~%")
  (loop for level from min-level upto max-level
	for value = (fe-value (solve-laplace problem level order :solver solver)
			      position)
	for new-approx = (vref (aref value 0) 0)
	and old-approx = nil then new-approx do
	(format t "~5D :  ~12,6,2E" level new-approx)
	(when old-approx (format t "   :  ~12,6,2E" (- new-approx old-approx)))
	(terpri)))

(defun check-p-convergence (problem min-order max-order
			    &key level position isopar (solver *lu-solver*))
  (format t "~%Dimension = ~D, level = ~D, midpoint = ~A:~%"
	  (dimension (domain problem)) level position)
  (format t "Order :  u_h(midpoint)  :  (u_h-u_2h)(midpoint)~%")
  (loop for order from min-order upto max-order
	for value =
	(fe-value (solve-laplace problem level order :solver solver :parametric
				 (and isopar (lagrange-mapping order)))
		  position)
	for new-approx = (vref (aref value 0) 0)
	and old-approx = nil then new-approx do
	(format t "~5D :  ~12,6,2E" order new-approx)
	(when old-approx (format t "   :  ~12,6,2E" (- new-approx old-approx)))
	(terpri)))


