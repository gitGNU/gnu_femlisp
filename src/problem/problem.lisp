;;; -*- mode: lisp; -*-

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; problem.lisp
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

(in-package :problem)

(defclass <problem> ()
  ()
  (:documentation "Base class for all problems."))

(defgeneric linear-p (problem)
  (:documentation "Predicate determining if a problem is linear or nonlinear."))

(defgeneric ensure-residual (problem blackboard)
  (:documentation "Ensures that the field :RESIDUAL is computed and that
the flag :RESIDUAL-P is set on the blackboard."))

(defgeneric ensure-solution (problem blackboard)
  (:documentation "Ensures that the field :SOLUTION is set on the
blackboard."))

(defmethod ensure-solution :around (problem blackboard)
  "This :around method takes the field :SOLUTION from the blackboard."
  (or (getbb blackboard :solution)
      (call-next-method)))

(defmethod ensure-solution (problem blackboard)
  (error "No initial guess for SOLUTION."))

(defmethod ensure-residual (problem blackboard)
  "This default method handles nonlinear problems by linearizing them and
computing the residual for it.  The resulting problem is additionally
stored in a field :LINEAR-PROBLEM."
  (when (linear-p problem)
    (error "No method ENSURE-RESIDUAL provided for this linear problem."))
  ;; linearize and compute residual
  (with-items (&key linearization solution) blackboard
    (setq linearization (linearize problem solution))
    (ensure-residual linearization blackboard)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; General problem solving
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass <solver> ()
  ((output :reader output :initform nil :initarg :output))
  (:documentation "The base class of linear, nonlinear and whatever
iterative solvers."))

(defgeneric solve (solver &optional blackboard)
  (:documentation "Solve a problem specified on the blackboard.  Returns a
modified blackboard.  The returned blackboard is guaranteed to contain at
least the fields :solution and :status.  :status is one of the values
:success or :failure.

SOLVE can also be called as (SOLVE blackboard) and will then try to figure
out a suitable solver itself."))

(defgeneric select-solver (problem blackboard)
  (:documentation "Choose a solver for PROBLEM.  The choice can also be
influenced by other data on the BLACKBOARD."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Linear and nonlinear problems
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Linear problems

(defclass <lse> (<problem>)
  ((matrix :reader matrix :initarg :matrix :initform nil)
   (rhs :reader rhs :initarg :rhs :initform nil))
  (:documentation "Standard form of a linear system of equations."))

(defun lse (&rest args)
  "Constructs a standard LSE."
  (apply #'make-instance '<lse> args))

(defmethod linear-p ((lse <lse>)) t)

(defgeneric linearize (problem solution)
  (:documentation "Linearize the nonlinear problem PROBLEM at the point
SOLUTION.  The result should be a linear problem."))

(defmethod linearize (problem solution)
  "This default method throws an error for nonlinear problems and is the
identity on linear problems."
  (if (linear-p problem)
      problem
      (error "No method LINEARIZE provided for this problem.")))

(defmethod ensure-solution ((lse <lse>) blackboard)
  (setf (getbb blackboard :solution)
	(make-domain-vector-for
	 (matrix lse) (multiplicity (rhs lse)))))

(defmethod ensure-residual ((lse <lse>) blackboard)
  (with-items (&key solution residual residual-p) blackboard
    (unless residual
      (setq residual (make-image-vector-for
		      (matrix lse) (multiplicity (rhs lse)))))
    (unless residual-p
      (copy! (rhs lse) residual)
      (gemm! -1.0 (matrix lse) solution 1.0 residual)
      (setq residual-p t))))

;;; Nonlinear problems

(defclass <nlse> (<problem>)
  ((linearization :initarg :linearization
    :documentation "A function linearizing the problem."))
  (:documentation "Class for nonlinear system of equations.  The
linearization contains a function returning a linear problem."))

(defmethod linear-p ((nlse <nlse>)) nil)

(defmethod linearize ((problem <nlse>) solution)
  (funcall (slot-value problem 'linearization) solution))

(defun nlse (&rest args)
  "Constructs a standard NLSE."
  (apply #'make-instance '<nlse> args))

;;; Testing: (problem::test-problem)

(defun test-problem ()
  (describe (lse :matrix #m(1.0) :rhs #m(1.0)))
  ;; problem for solving x^2=2 by Newton's method
  (let ((nlse
	 (nlse :linearization
	       #'(lambda (solution)
		   (let ((x (mref solution 0 0)))
		     (lse :matrix (make-real-matrix (vector (* 2.0 x)))
			  :rhs (make-real-matrix (vector (+ 2.0 (* x x)))))))))
	(bb (blackboard :solution #m(1.4))))
    (ensure-residual nlse bb)
    bb)
  )
(fl.tests:adjoin-test 'test-problem)