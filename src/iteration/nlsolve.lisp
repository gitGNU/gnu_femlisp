;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; nlsolve.lisp - Provide nonlinear solvers
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

(in-package :fl.iteration)

(file-documentation
 "This file provides some solvers for nonlinear problems.")

(defclass <nonlinear-solver> (<discrete-iterative-solver>)
  ((linear-solver
    :reader linear-solver :initarg :linear-solver
    :documentation "The linear solver for solving the linearization."))
  (:documentation "Class for general nonlinear iterative solvers."))

(defmethod inner-iteration ((nlsolve <nonlinear-solver>))
  (and (slot-boundp nlsolve 'linear-solver)
       (linear-solver nlsolve)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Newton iteration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <newton> (<nonlinear-solver>)
  ()
  (:documentation "Class for the Newton iteration."))

;;; The following does not work, because blackboard is not available.  This
;;; might suggest a change in interface such that all solvers and
;;; strategies should carry their own blackboards around as a slot.
#+(or)
(defmethod slot-unbound (class (newton <newton>) slot)
  (case slot
    (linear-solver (select-linear-solver newton blackboard))))

(defmethod next-step ((newton <newton>) blackboard)
  "Simply calls the linear solver on the linearized problem."
  (with-items (&key solution residual residual-p linearization)
      blackboard
    (unless (slot-boundp newton 'linear-solver)
      (setf (slot-value newton 'linear-solver)
	    (select-linear-solver linearization blackboard)))
    (solve (linear-solver newton)
	   (blackboard :problem linearization :solution solution
		       :residual residual :residual-p t))
    ;; and since the linear residual is not correct for the nonlinear
    ;; problem:
    (setq residual-p nil)))

(defmethod select-solver ((problem <nonlinear-problem>) blackboard)
  (make-instance
   '<newton>
   :success-if `(and (> :step 1) (> :step-reduction 0.5))
   :failure-if '(and (> :step 1) (> :step-reduction 1.0) (> :defnorm 1.0e-5))
   :output (getbb blackboard :output)))


;;;; Testing

(defun test-nlsolve ()
  (let ((nlse
	 (nlse :linearization
	       #'(lambda (solution)
		   (let ((x (mref solution 0 0)))
		     (lse :matrix (make-real-matrix (vector (* 2.0 x)))
			  :rhs (make-real-matrix (vector (+ 2.0 (* x x))))))))))
    (solve (make-instance
	    '<newton> :output t
	    :success-if '(and (> :step 2) (> :step-reduction 0.9))
	    :observe (append *discrete-iterative-solver-observe*
			     `(("                 solution" "~25,15,2E"
				,#'(lambda (blackboard)
				    (with-items (&key solution) blackboard
				      (vref solution 0)))))))
	   (blackboard :problem nlse :solution #m(1.0))))
  )
  
;;;; Testing: (test-nlsolve)
(fl.tests:adjoin-test 'test-nlsolve)