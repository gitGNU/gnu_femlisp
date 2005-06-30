;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; solve.lisp - Provide iterative solvers
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; General solvers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <solver> ()
  ((output :reader output :initarg :output))
  (:documentation "The base class of linear, nonlinear and whatever
iterative solvers."))

(defgeneric solve (solver &optional blackboard)
  (:documentation "Solve a problem specified on the blackboard.  Returns a
modified blackboard.  The returned blackboard is guaranteed to contain at
least the fields :solution and :status.  :status is one of the values
:success or :failure.

SOLVE can also be called as (SOLVE blackboard) and will then try to figure
out a suitable solver itself."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Iterative solvers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <iterative-solver> (<iteration> <solver>)
  ()
  (:documentation "Base class of all iterative solvers and solution
strategies."))

(defmethod solve ((itsol <iterative-solver>) &optional blackboard)
  (dbg :iter "Iterating: ~S" itsol)
  (iterate itsol blackboard))

;;; Problem-dependent solving

(defgeneric select-solver (object blackboard)
  (:documentation "Selects a solver for OBJECT.  OBJECT is usually a
problem with certain characteristics."))

(defgeneric select-linear-solver (object blackboard)
  (:documentation "Selects a linear solver for OBJECT.  OBJECT is usually a
matrix or a linear problem with certain characteristics."))

(defmethod select-solver :around (object blackboard)
  "If a solver is on the blackboard, use it."
  (declare (ignore object))
  (or (getbb blackboard :solver) (call-next-method)))

(defmethod select-solver ((problem <problem>) blackboard)
  "This method does a more specific search for linear problems by calling
@function{select-linear-solver}."
  (if (linear-p problem)
      (select-linear-solver problem blackboard)
      (error "No solver for this problem known.")))

(defvar *select-linear-solver* nil
  "If this variable is non-nil, it is funcalled for selecting a linear
solver.  An alternative way of handling this could be given by Pascal
Costanza's @code{AspectL}.")

(defmethod select-linear-solver :around (object blackboard)
  "If a linear solver is on the blackboard, use it."
  (or (and *select-linear-solver*
	   (funcall *select-linear-solver* object blackboard))
      (getbb blackboard :solver) (call-next-method)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Iterative solvers for discrete problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *discrete-iterative-solver-observe*
  (list (list "Step" "~4D"
	      #'(lambda (blackboard) (getbb blackboard :step)))
	(list "||residual||" "~12,5,2E"
	      #'(lambda (blackboard) (getbb blackboard :defnorm)))
	(list "res[k]/res[k-1]" "   ~12,5,2E"
	      #'(lambda (blackboard)
		  (getbb blackboard :step-reduction))))
  "Standard observe quantities for iterative solvers.")

(defclass <discrete-iterative-solver> (<iterative-solver>)
  ((residual-norm :reader residual-norm :initform #'norm
		  :initarg :residual-norm)
   (observe :initform *discrete-iterative-solver-observe*))
  (:documentation "The base class of solvers for discrete linear and
nonlinear problems."))

(defmethod initially ((itsol <discrete-iterative-solver>) blackboard)
  "Ensure an initial guess for the solution after other initializations
have been done."
  (call-next-method)
  (with-items (&key problem defnorm initial-defnorm residual-p) blackboard
    (ensure-solution problem blackboard)
    (setq defnorm nil initial-defnorm nil residual-p nil)))

(defun safe-divide-by-zero (a b)
  (if (zerop a)
      (if (zerop b) :undefined 0.0)
      (if (zerop b) :infinity (/ a b))))

(defmethod intermediate :before ((itsolve <discrete-iterative-solver>) blackboard)
  "Before printing information in the main method we ensure that the defect
and its norm are up-to-date."
  (with-items (&key problem residual step defnorm previous-defnorm initial-defnorm
		    reduction step-reduction)
      blackboard
    (ensure-residual problem blackboard)
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

(defmethod finally ((itsolve <discrete-iterative-solver>) blackboard)
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
  (make-instance '<discrete-iterative-solver>)
  )
(fl.tests:adjoin-test 'test-solve)
