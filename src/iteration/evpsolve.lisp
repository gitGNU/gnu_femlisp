;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; evpsolve.lisp - Provide solvers for eigenvalue problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Copyright (C) 2004 Nicolas Neuss, University of Heidelberg.
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

(in-package :fl.iteration)

(file-documentation
 "This file provides solvers for eigenvalue problems.")

(defclass <evp-solver> (<nonlinear-solver>)
  ()
  (:documentation "General class for a solver of eigenvalue problems."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Newton iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <wielandt-iteration> (<evp-solver>)
  ()
  (:documentation "Wielandt iteration."))

(defmethod next-step ((wielandt <wielandt-iteration>) blackboard)
  "Simply calls the linear solver on the linearized problem."
  (with-items (&key problem solution residual residual-p
		    linearization linear-solver-failure step)
      blackboard
    (unless (slot-boundp wielandt 'linear-solver)
      (setf (slot-value wielandt 'linear-solver)
	    (?1 (lu-solver)
		(select-linear-solver linearization blackboard))))
    (copy! solution (rhs linearization))
    (handler-case
       (solve (linear-solver wielandt)
	      (blackboard :problem linearization :solution solution
			  :residual residual :residual-p residual-p))
     (arithmetic-error ()
       ;; In the case of an arithmetic error we have found an eigenvalue.
       ;; In practice, this will rarely occur.  Note that if this occurs in
       ;; step 1, the eigenvector may be completely wrong.
       (setq linear-solver-failure t)))
    (let ((mass (mass problem solution))
	  (energy (energy problem solution)))
      (setf (unbox (slot-value problem 'lambda)) (/ energy mass))
      (scal! (/ (sqrt mass)) solution)
      (dbg :iter "mass=~A energy=~A lambda=~A x=~A~%"
	   mass energy (unbox (slot-value problem 'lambda)) solution))
    (setq residual-p nil)))

(defmethod select-solver ((problem <evp>) blackboard)
  (make-instance
   '<wielandt-iteration>
   :success-if (or (getbb blackboard :success-if)
		   '(and (> :step 1) (> :step-reduction 0.5)))
   :failure-if (or (getbb blackboard :failure-if)
		   '(and (> :step 1) (> :step-reduction 1.0) (> :defnorm 1.0e-5)))))

;;;; Testing

(defun test-evpsolve ()
  (let* ((dim 3)
	 (evp (make-instance
	       '<ls-evp> :A (laplace-full-matrix dim) :B (eye dim)))
	 (bb (blackboard :problem evp :output t)))
    (solve bb)
    (unbox (slot-value evp 'lambda)))
  )

;;;; Testing: (test-evpsolve)
(fl.tests:adjoin-test 'test-evpsolve)
