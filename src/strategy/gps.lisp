;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; gps.lisp - General problem solving strategy 
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

(in-package :strategy)

;;;; Problem analysis and solution with an automatic choice of solver.
;;;; Very first version.

(defmethod solve ((blackboard list) &optional dummy)
  "Implements the simplest interface for problem solving.  Only a
blackboard describing the problem is passed and the solution method is
open."
  (when (or dummy (not (eq (car blackboard) :blackboard)))
    (error "Syntax: (SOLVE blackboard) or (SOLVE strategy blackboard)."))
  (with-items (&key strategy problem matrix rhs fe-class solver
		    success-if output plot-mesh)
      blackboard
      (cond
	(strategy			; there is a strategy on the blackboard
	 (solve strategy blackboard))
	(problem			; there is a PDE problem on the blackboard
	 (let* ((dim (dimension (domain problem)))
		(nr-comps (typecase problem
			    (cdr::<cdr-problem> 1)
			    (elasticity::<elasticity-problem> dim)
			    (navier-stokes::<navier-stokes-problem> (1+ dim))))
		(order (if (= nr-comps 1)
			   (if (<= dim 2) 6 4)
			   (if (= nr-comps dim)
			       (if (<= dim 2) 5 4)
			       (if (<= dim 2) 4 3))))
		(disc (typecase problem
			(cdr::<cdr-problem> (lagrange-fe order))
			(elasticity::<elasticity-problem> (lagrange-fe order :nr-comps dim))
			(navier-stokes::<navier-stokes-problem>
			 (navier-stokes-fe::navier-stokes-lagrange-fe order dim 1)))))
	   (let ((strategy
		  (make-instance
		   '<stationary-fe-strategy>
		   :fe-class (or fe-class disc) :solver (or solver *lu-solver*)
		   :indicator (make-instance '<uniform-refinement-indicator>)
		   :success-if success-if :output output
		   :plot-mesh plot-mesh)))
	     (solve strategy blackboard))))
	((and matrix rhs)		; there is a linear problem on the blackboard
	 (solve (or solver *lu-solver*) blackboard))
	(t (error "SOLVE does not know about this problem.")))))


