;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; solve.lisp
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

(in-package :iteration)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; General solver class
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <solver> ()
  ((output :reader output :initform nil :initarg :output))
  (:documentation "The base class of linear, nonlinear and whatever
iterative solvers."))

(defgeneric solve (solver blackboard)
  (:documentation "Solve a problem specified on the blackboard.  Returns a
modified blackboard.  The returned blackboard is guaranteed to contain at
least the fields :solution and :status.  :status may have the values
:success or :failure."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Iterative solvers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *iterative-solver-observe*
  (list (list "Step" "~4D"
	      #'(lambda (blackboard) (getbb blackboard :step)))
	(list "||residual||" "~12,5,2E"
	      #'(lambda (blackboard) (getbb blackboard :defnorm)))
	(list "res[k]/res[k-1]" "   ~12,5,2E"
	      #'(lambda (blackboard)
		  (getbb blackboard :step-reduction))))
  "Standard observe quantities for iterative solvers.")

(defclass <iterative-solver> (<iteration> <solver>)
  ((iteration :reader iteration :initarg :iteration
	      :documentation "The inner iteration.  This iteration should
perform reasonably efficient.")
   (residual-norm :reader residual-norm :initform #'norm
		  :initarg :residual-norm)
   (observe :initform *iterative-solver-observe*))
  (:documentation "The base class of linear, nonlinear and whatever
iterative solvers."))

(defgeneric initial-guess (itsol blackboard)
  (:documentation "Generates an initial guess for a solution to the
problem.  Pre-condition: problem.  Post-condition: solution."))

(defgeneric ensure-residual (itsol blackboard)
  (:documentation "Ensures that the residual is in a valid state."))

(defmethod name ((itsolve <iterative-solver>))
  "The name includes also information on the inner iteration."
  (format nil "~A (~A)"
	  (class-name (class-of itsolve))
	  (class-name (class-of (iteration itsolve)))))

(defmethod initially ((itsol <iterative-solver>) blackboard)
  "Ensure an initial guess after other initializations have been done."
  (dbg :iter "Initially: blackboard = ~A" blackboard)
  (call-next-method)
  (initial-guess itsol blackboard)
  (with-items (&key defnorm initial-defnorm residual-p) blackboard
    (setq defnorm nil initial-defnorm nil residual-p nil)))

(defmethod solve ((itsol <iterative-solver>) blackboard)
  (iterate itsol blackboard))

(defun safe-divide-by-zero (a b)
  (if (zerop a)
      (if (zerop b) :undefined 0.0d0)
      (if (zerop b) :infinity (/ a b))))

(defmethod intermediate :before ((itsolve <iterative-solver>) blackboard)
  "Before printing information in the main method we ensure that the defect
and its norm are up-to-date."
  (with-items (&key solution matrix rhs residual
		    step defnorm previous-defnorm initial-defnorm
		    reduction step-reduction)
      blackboard
    (ensure-residual itsolve blackboard)
    ;; set new residual norm
    (setq previous-defnorm defnorm)
    (setq defnorm (funcall (residual-norm itsolve) residual))
    (unless initial-defnorm (setq initial-defnorm defnorm))
    ;; compute derived quantities
    (setq reduction (safe-divide-by-zero defnorm initial-defnorm))
    (setq step-reduction
	  (if (zerop step)
	      :undefined
	      (safe-divide-by-zero defnorm previous-defnorm)))))

(defmethod next-step ((itsolve <iterative-solver>) blackboard)
  "Stepping for a linear solver."
  (with-items (&key solution matrix rhs residual residual-p step iterator)
      blackboard
    (unless iterator ; initialize the inner iteration if necessary
      (dbg :iter "Linear solving: making iterator")
      (setq iterator (make-iterator (iteration itsolve) matrix))
      (awhen (slot-value iterator 'initialize)
	(funcall it solution rhs residual)))
    (funcall (slot-value iterator 'iterate)
	     solution rhs residual)
    (setq residual-p (slot-value iterator 'residual-after))))

(defmethod finally ((itsolve <iterative-solver>) blackboard)
  "Put convergence rate into report."
  (call-next-method)
  (with-items (&key report step status reduction) blackboard
    (setf (getbb report :convergence-rate)
	  (if (and (plusp step) (numberp reduction))
	      (expt reduction (/ step))
	      :undefined)))
  (dbg :iter "Finally: blackboard = ~A" blackboard))



;;;; Testing: (test-solve)
(defun test-solve ()
  (make-instance '<iterative-solver>)
  )
(tests::adjoin-femlisp-test 'test-solve)
